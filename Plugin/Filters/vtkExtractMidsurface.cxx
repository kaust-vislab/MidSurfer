#include "vtkExtractMidsurface.h"
#include <vtkImageData.h>
#include <vtkGradientFilter.h>
#include <vtkArrayCalculator.h>
#include <vtkPassSelectedArrays.h>
#include <vtkDataArraySelection.h>
#include <vtkPassArrays.h>
#include <vtkAppendFilter.h>
#include <vtkExtractVectorComponents.h>
#include <vtkDataSetAttributes.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkVector.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkStreamTracer.h>
#include <vtkLineSource.h>
#include <vtkLine.h>
#include <vtkExtractVOI.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include <vtkContourFilter.h>
#include <vtkPointLocator.h>
#include <vtkNew.h>
#include <vtkFeatureEdges.h>
#include <vtkCenterOfMass.h>
#include <vtkGeometryFilter.h>
#include <vtkSmartPointer.h>

#include "vtkSignedDistanceField.h"

#include <chrono>
#include <numbers>

#define SHFT2(a, b, c) \
	(a) = (b);         \
	(b) = (c);
#define SHFT3(a, b, c, d) \
	(a) = (b);            \
	(b) = (c);            \
	(c) = (d);

vtkStandardNewMacro(vtkExtractMidsurface);

//----------------------------------------------------------------------------
vtkExtractMidsurface::vtkExtractMidsurface() : InputArray(nullptr)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::vtkExtractMidsurface() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::vtkExtractMidsurface() END");
}

vtkExtractMidsurface::~vtkExtractMidsurface()
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::~vtkExtractMidsurface() BEGIN");

	this->SetInputArray(nullptr);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::~vtkExtractMidsurface() END");
}

int vtkExtractMidsurface::FillOutputPortInformation(int, vtkInformation *info)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::FillOutputPortInformation() BEGIN");

	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::FillOutputPortInformation() END");

	return 1;
}

int vtkExtractMidsurface::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

	int extent[6];
	double origin[3];
	double spacing[3];

	inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);
	inInfo->Get(vtkDataObject::SPACING(), spacing);
	inInfo->Get(vtkDataObject::ORIGIN(), origin);

	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
	outInfo->Set(vtkDataObject::SPACING(), spacing, 3);
	outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);

	return 1;
}

int vtkExtractMidsurface::RequestData(vtkInformation *vtkNotUsed(request),
									  vtkInformationVector **inputVector,
									  vtkInformationVector *outputVector)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::RequestData() BEGIN");

	const auto timer_start{std::chrono::steady_clock::now()};

	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	auto input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	auto output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkNew<vtkImageData> image;
	image->DeepCopy(input); // for initilaization

	int dims[3];
	double spacing[3];
	double bounds[6];

	image->GetDimensions(dims);
	image->GetSpacing(spacing);
	image->GetBounds(bounds); // not used...
	
	this->SetProgressText("Computing step size.");
	this->UpdateProgress(0.000001);

	if (this->AutomaticStepSize)
		this->IntegrationStep = std::numbers::sqrt2 * spacing[2];

	this->SetProgressText("Computing smoothed signed distance map.");
	this->UpdateProgress(0.000002);

	if (this->SmoothInput == SMOOTH_INPUT_GAUSSIAN_3D)
	{
		ComputeGaussianSmoothing(image);
	}
	else if (this->SmoothInput == SMOOTH_INPUT_SDF_3D)
	{
		ComputeSmoothSignedDistanceMap(image);
	}

	vtkNew<vtkAppendPolyData> append;

	if (dims[2] == 1)
		vtkErrorMacro("ExtractMidsurface only works on a volume (dim z > 1).");
	else
	{
		ExtractMidsurface(image, append, dims);
	}

	output->ShallowCopy(append->GetOutput());

	const auto timer_end{std::chrono::steady_clock::now()};
	const std::chrono::duration<double> elapsed_time(timer_end - timer_start);
	vtkLog(INFO, "Elapsed time: " << std::to_string(elapsed_time.count()));
	// vtkLog(INFO, "Elapsed time: " << elapsed_time); // should work but doesn't

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::RequestData() END");

	return 1;
}

void vtkExtractMidsurface::ExtractMidsurface(vtkImageData *image, vtkAppendPolyData *append, int *dims)
{
	for (int z = 0; z < dims[2]; z++)
	{
		vtkNew<vtkExtractVOI> slice;
		slice->SetInputData(image);
		slice->SetVOI(0, dims[0] - 1, 0, dims[1] - 1, z, z);
		slice->Update();

		vtkNew<vtkExtractCenterLine> centerline;
		centerline->SetInputData(slice->GetOutput());
		centerline->SetInputArray(this->InputArray);
		centerline->SetAutomaticStepSize(false);
		centerline->SetDistanceType(this->DistanceType);
		centerline->SetFieldType(this->FieldType);
		centerline->SetGoldenSectionSearch(this->GoldenSectionSearch);
		centerline->SetIntegrationStep(this->IntegrationStep);
		centerline->SetRadiusFactors(this->RadiusFactors);
		centerline->SetRadiusFactorsSDF(this->RadiusFactorsSDF);
		centerline->SetResultType(this->ResultType);
		centerline->SetShapeDetection(this->ShapeDetection);
		centerline->SetSmoothInput((this->SmoothInput == ESmoothingMethod::SMOOTH_INPUT_GAUSSIAN_3D || this->SmoothInput == ESmoothingMethod::SMOOTH_INPUT_SDF_3D) ? ESmoothingMethod::SMOOTH_INPUT_OFF : this->SmoothInput);
		centerline->SetStandardDeviations(this->StandardDeviations);
		centerline->SetStandardDeviationsSDF(this->StandardDeviationsSDF);
		centerline->SetThreshold(this->Threshold);
		centerline->SetTolerance(this->Tolerance);
		centerline->SetConnectivity(this->Connectivity);
		centerline->Update();

		append->AddInputConnection(centerline->GetOutputPort());

		this->SetProgressText(("Processing slice " + std::to_string(z) + "/" + std::to_string(dims[2])).c_str());
		this->UpdateProgress(static_cast<double>(z) / static_cast<double>(dims[2]));
	}

	append->Update();
}

void vtkExtractMidsurface::ComputeSmoothSignedDistanceMap(vtkImageData *image)
{
	vtkLog(INFO, "Smoothing input with SDF");

	vtkNew<vtkSignedDistanceField> sdf;
	sdf->SetInputData(image);
	sdf->SetInputArray(this->InputArray);
	sdf->SetDistanceType(this->DistanceType);
	sdf->SetFieldType(this->FieldType);
	sdf->SetSmoothing(this->Smoothing);
	sdf->SetManualSmoothing(this->ManualSmoothing);
	sdf->SetRadiusFactors(this->SDFRadiusFactors);
	sdf->SetStandardDeviations(this->SDFStandardDeviations);
	sdf->SetSigmaDivisor(this->SigmaDivisor);
	sdf->SetTestingRadius(this->TestingRadius);
	sdf->SetNormalize(this->Normalize);
	sdf->SetFinalResultName("smooth");
	sdf->Update();

	image->GetPointData()->AddArray(sdf->GetOutput()->GetPointData()->GetArray("smooth"));
}

void vtkExtractMidsurface::ComputeGaussianSmoothing(vtkImageData *image)
{
	// copy the input array to double array
	vtkNew<vtkArrayCalculator> calc;
	calc->SetInputData(image);
	calc->AddScalarArrayName(this->InputArray);
	calc->SetFunction(this->InputArray);
	calc->SetResultArrayName("smooth");
	calc->SetResultArrayType(VTK_DOUBLE);
	calc->Update();

	vtkNew<vtkImageGaussianSmooth> smooth;
	smooth->SetInputData(calc->GetOutput());
	smooth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "smooth");
	smooth->SetRadiusFactors(this->RadiusFactors);
	smooth->SetStandardDeviations(this->StandardDeviations);
	smooth->Update();

	image->GetPointData()->AddArray(smooth->GetOutput()->GetPointData()->GetArray("smooth"));
}
