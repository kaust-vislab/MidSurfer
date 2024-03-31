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
#include "vtkZipperTriangulation.h"

#include <chrono>
#include <numbers>
#include <algorithm>

#define SHFT2(a, b, c) \
	(a) = (b);         \
	(b) = (c);
#define SHFT3(a, b, c, d) \
	(a) = (b);            \
	(b) = (c);            \
	(c) = (d);

vtkStandardNewMacro(vtkExtractMidsurface);

//----------------------------------------------------------------------------
vtkExtractMidsurface::vtkExtractMidsurface() : InputArray(nullptr), NumberOfLabels(0)
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

	double range[2];
	input->GetPointData()->GetArray(this->InputArray)->GetRange(range);
	this->NumberOfLabels = static_cast<int>(range[1]);
	vtkLog(INFO, "Number of labels: " << this->NumberOfLabels);

	vtkNew<vtkAppendPolyData> append;

	for (int labelId = 1; labelId <= this->NumberOfLabels; labelId++)
	{
		vtkNew<vtkArrayCalculator> calc;
		calc->SetInputData(input);
		calc->AddScalarArrayName(this->InputArray);
		calc->SetFunction((std::string(this->InputArray) + " == " + std::to_string(labelId)).c_str());
		calc->SetResultArrayName("input_label_mask");
		calc->SetResultArrayType(VTK_DOUBLE);
		calc->Update();

		vtkImageData *image = vtkImageData::SafeDownCast(calc->GetOutput());

		int dims[3];
		int extent[6];
		double spacing[3];
		double bounds[6];

		image->GetDimensions(dims);
		image->GetSpacing(spacing);
		image->GetBounds(bounds); // not used...
		image->GetExtent(extent);

		this->SetProgressText(("Computing VOI (" + std::to_string(labelId) + "/" + std::to_string(this->NumberOfLabels) + ").").c_str());
		this->UpdateProgress(static_cast<double>(labelId - 1) / static_cast<double>(this->NumberOfLabels));

		int labelExtent[6];
		FindLabelExtent(labelExtent, extent, image);

		vtkNew<vtkExtractVOI> extractVOI;
		extractVOI->SetInputData(image);
		extractVOI->SetVOI(labelExtent[0], labelExtent[1], labelExtent[2], labelExtent[3], labelExtent[4], labelExtent[5]);
		extractVOI->Update();

		auto voi = extractVOI->GetOutput();

		if (this->AutomaticStepSize)
			this->IntegrationStep = std::numbers::sqrt2 * spacing[2];

		if (this->SmoothInput == SMOOTH_INPUT_GAUSSIAN_3D)
		{
			ComputeGaussianSmoothing(voi);
		}
		else if (this->SmoothInput == SMOOTH_INPUT_SDF_3D)
		{
			ComputeSmoothSignedDistanceMap(voi);
		}

		if (dims[2] == 1)
			vtkErrorMacro("ExtractMidsurface only works on a volume (dim z > 1).");
		else
		{
			vtkNew<vtkPolyData> mesh;
			ExtractMidsurface(voi, mesh, labelId);
			append->AddInputData(mesh);
		}
	}

	append->Update();
	output->ShallowCopy(append->GetOutput());

	const auto timer_end{std::chrono::steady_clock::now()};
	const std::chrono::duration<double> elapsed_time(timer_end - timer_start);
	vtkLog(INFO, "Elapsed time: " << std::to_string(elapsed_time.count()));
	// vtkLog(INFO, "Elapsed time: " << elapsed_time); // should work but doesn't

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface::RequestData() END");

	return 1;
}

void vtkExtractMidsurface::FindLabelExtent(int *labelExtent, int *extent, vtkImageData *image)
{
	labelExtent[0] = VTK_INT_MAX;
	labelExtent[1] = VTK_INT_MIN;
	labelExtent[2] = VTK_INT_MAX;
	labelExtent[3] = VTK_INT_MIN;
	labelExtent[4] = VTK_INT_MAX;
	labelExtent[5] = VTK_INT_MIN;

	for (int z = extent[4]; z < extent[5]; z++)
		for (int y = extent[2]; y < extent[3]; y++)
			for (int x = extent[0]; x < extent[1]; x++)
			{
				if (*((double *)(image->GetScalarPointer(x, y, z))) > 0)
				{
					if (x < labelExtent[0])
						labelExtent[0] = x;
					if (x > labelExtent[1])
						labelExtent[1] = x;
					if (y < labelExtent[2])
						labelExtent[2] = y;
					if (y > labelExtent[3])
						labelExtent[3] = y;
					if (z < labelExtent[4])
						labelExtent[4] = z;
					if (z > labelExtent[5])
						labelExtent[5] = z;
				}
			}

	labelExtent[0] = std::max({extent[0], labelExtent[0] - this->LabelExtentBorder});
	labelExtent[1] = std::min({extent[1], labelExtent[1] + this->LabelExtentBorder});
	labelExtent[2] = std::max({extent[2], labelExtent[2] - this->LabelExtentBorder});
	labelExtent[3] = std::min({extent[3], labelExtent[3] + this->LabelExtentBorder});
	labelExtent[4] = std::max({extent[4], labelExtent[4] - this->LabelExtentBorder});
	labelExtent[5] = std::min({extent[5], labelExtent[5] + this->LabelExtentBorder});

	vtkLog(INFO, "VOI: [" << labelExtent[0] << ", " << labelExtent[1] << ", "
						  << labelExtent[2] << ", " << labelExtent[3] << ", "
						  << labelExtent[4] << ", " << labelExtent[5] << "]");
}

void vtkExtractMidsurface::ExtractMidsurface(vtkImageData *image, vtkPolyData *mesh, int labelId)
{
	int extent[6];
	image->GetExtent(extent);

	vtkNew<vtkAppendPolyData> append;

	for (int z = extent[4]; z < extent[5]; z++)
	{
		vtkNew<vtkExtractVOI> slice;
		slice->SetInputData(image);
		// slice->SetVOI(0, dims[0] - 1, 0, dims[1] - 1, z, z);
		slice->SetVOI(extent[0], extent[1], extent[2], extent[3], z, z);
		slice->Update();

		vtkNew<vtkExtractCenterLine> centerline;
		centerline->SetInputData(slice->GetOutput());
		centerline->SetInputArray("input_label_mask");
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

		this->SetProgressText(("Processing slice " + std::to_string(z - extent[4]) + "/" + std::to_string(extent[5] - extent[4])).c_str());
		this->UpdateProgress(static_cast<double>(labelId - 1) / static_cast<double>(this->NumberOfLabels) +
							 static_cast<double>(z - extent[4]) / (static_cast<double>(extent[5] - extent[4]) * static_cast<double>(this->NumberOfLabels)));
	}

	append->Update();

	auto num = append->GetOutput()->GetNumberOfPoints();

	vtkNew<vtkDoubleArray> surfaceId;
	surfaceId->SetName("SurfaceId");
	surfaceId->SetNumberOfTuples(num);
	for (unsigned int i = 0; i < num; i++)
	{
		surfaceId->InsertTuple1(i, labelId);
	}

	append->GetOutput()->GetPointData()->AddArray(surfaceId);

	if (this->ResultType == Midsurfacer::Tools::RESULT_TYPE_TRIANGULATION)
	{
		vtkNew<vtkZipperTriangulation> zipper;
		zipper->SetInputConnection(append->GetOutputPort());
		zipper->SetZipperAlpha(this->ZipperAlpha);
		zipper->Update();
		mesh->DeepCopy(zipper->GetOutput());
	}
	else
		mesh->DeepCopy(append->GetOutput());
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
