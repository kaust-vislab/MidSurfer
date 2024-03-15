#include "vtkExtractMidsurface3D.h"
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
#include <vtkDataArray.h>
#include <vtkVector.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkStreamTracer.h>
#include <vtkLineSource.h>
#include <vtkLine.h>
#include <vtkExtractVOI.h>
#include <vtkThreshold.h>
#include <vtkConnectivityFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageCast.h>

vtkStandardNewMacro(vtkExtractMidsurface3D);

//----------------------------------------------------------------------------
vtkExtractMidsurface3D::vtkExtractMidsurface3D() : InputArray(nullptr)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractMidsurface3D::vtkExtractMidsurface3D() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->a = create_matrix<double>(3, 3);
	this->w = create_matrix<double>(3, 3);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), ":vtkExtractMidsurface3D:vtkExtractMidsurface3D() END");
}

vtkExtractMidsurface3D::~vtkExtractMidsurface3D()
{
	free_matrix(this->w);
	free_matrix(this->a);

	this->SetInputArray(nullptr);
}

int vtkExtractMidsurface3D::FillOutputPortInformation(int, vtkInformation *info)
{
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
	return 1;
}

template <class TReal>
TReal **vtkExtractMidsurface3D::create_matrix(long nrow, long ncol)
{
	typedef TReal *TRealPointer;
	TReal **m = new TRealPointer[nrow];

	TReal *block = (TReal *)calloc(nrow * ncol, sizeof(TReal));
	m[0] = block;
	for (int row = 1; row < nrow; ++row)
	{
		m[row] = &block[row * ncol];
	}
	return m;
}

/* free a TReal matrix allocated with matrix() */
template <class TReal>
void vtkExtractMidsurface3D::free_matrix(TReal **m)
{
	free(m[0]);
	delete[] m;
}

void vtkExtractMidsurface3D::ExtractMidsurface(vtkImageData *input, vtkAppendPolyData *append)
{
	this->StartPoints.clear();

	vtkNew<vtkImageCast> castFilter;
	castFilter->SetOutputScalarTypeToDouble();
	castFilter->SetInputData(input);
	castFilter->Update();

	vtkNew<vtkImageGaussianSmooth> smooth;
	smooth->SetInputData(castFilter->GetOutput());
	smooth->SetRadiusFactors(this->RadiusFactors);
	smooth->SetStandardDeviations(this->StandardDeviations);
	smooth->Update();

	vtkNew<vtkThreshold> thresh;
	thresh->SetInputConnection(smooth->GetOutputPort());
	thresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, this->InputArray);
	thresh->SetLowerThreshold(0.5);
	thresh->SetUpperThreshold(10.0);
	thresh->Update();

	vtkNew<vtkConnectivityFilter> conn;
	conn->SetInputConnection(thresh->GetOutputPort());
	conn->SetExtractionModeToAllRegions();
	conn->ColorRegionsOn();
	conn->Update();

	vtkLog(INFO, "found " << conn->GetNumberOfExtractedRegions() << " connected regions");

	for (int k = 0; k < conn->GetNumberOfExtractedRegions(); k++)
	{
		vtkLog(INFO, "Processing region " << k);

		vtkNew<vtkThreshold> threshReg;
		threshReg->SetInputConnection(conn->GetOutputPort());
		threshReg->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
		threshReg->SetLowerThreshold(static_cast<double>(k));
		threshReg->SetUpperThreshold(static_cast<double>(k));
		threshReg->Update();

		auto region = threshReg->GetOutput();

		vtkIdType idMax = FindMaximumLocation(region->GetPointData()->GetArray("smooth"));
		if (idMax == -1)
		{
			vtkLog(INFO, "No maximum found.");
		}
		else
		{
			double *max = region->GetPoint(idMax);

			std::array<double, 3> arr{{max[0], max[1], max[2]}};
			this->StartPoints.push_back(arr);

			for (double d = 1.0; d < 5.5; d += 1.0)
			{
				std::array<double, 3> arr1{{max[0], max[1], max[2] - d}};
				std::array<double, 3> arr2{{max[0], max[1], max[2] + d}};
				this->StartPoints.push_back(arr1);
				this->StartPoints.push_back(arr2);
			}

			vtkLog(INFO, "Inserting max: " << max[0] << ", " << max[1] << ", " << max[2]);
		}
	}

	// Compute gradients
	auto grad = vtkSmartPointer<vtkGradientFilter>::New();
	grad->SetInputConnection(smooth->GetOutputPort());
	grad->SetResultArrayName("Gradient");
	grad->Update();

	// Compute tensors
	auto tens = vtkSmartPointer<vtkGradientFilter>::New();
	tens->SetInputConnection(grad->GetOutputPort());
	tens->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, "Gradient");
	tens->SetResultArrayName("T");
	tens->Update();

	auto tensArr = tens->GetOutput()->GetPointData()->GetArray("T");
	auto val = smooth->GetOutput()->GetPointData()->GetArray("ImageFile");

	for (std::array<double, 3> arr : this->StartPoints)
	{
		double p[3];
		vtkMath::Assign(arr, p);

		vtkLog(INFO, "Extracting line from " << p[0] << ", " << p[1] << ", " << p[2]);

		vtkNew<vtkPoints> points;
		vtkNew<vtkCellArray> lines;

		points->InsertNextPoint(p);

		int k = AppendPoints(p, 1, 1.0, input, points, lines, tensArr, val); // forward
		vtkMath::Assign(arr, p);
		AppendPoints(p, k, -1.0, input, points, lines, tensArr, val); // backward

		vtkNew<vtkPolyData> centerLine;
		centerLine->SetPoints(points);
		centerLine->SetLines(lines);

		append->AddInputData(centerLine);
	}
}

int vtkExtractMidsurface3D::RequestData(vtkInformation *vtkNotUsed(request),
									  vtkInformationVector **inputVector,
									  vtkInformationVector *outputVector)
{
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	auto input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	auto output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	int dims[3];
	double spacing[3];
	double bounds[6];

	input->GetDimensions(dims);
	input->GetSpacing(spacing);
	input->GetBounds(bounds);

	vtkNew<vtkAppendPolyData> append;

	ExtractMidsurface(input, append);

	append->Update();
	output->ShallowCopy(append->GetOutput());

	return 1;
}

vtkIdType vtkExtractMidsurface3D::FindMaximumLocation(vtkDataArray *arr)
{
	if (arr->GetNumberOfTuples() == 0)
	{
		vtkLog(INFO, "No array found.");
		return -1;
	}

	double val = arr->GetTuple1(0);
	vtkIdType idMax = 0;

	for (vtkIdType id = 0; id < arr->GetNumberOfTuples(); id++)
	{
		if (val < arr->GetTuple1(id))
		{
			val = arr->GetTuple1(id);
			idMax = id;
		}
	}

	return idMax;
}

double vtkExtractMidsurface3D::InsertNextPoint(double *p, double *v, double direction, double *vold, vtkPoints *points)
{

	if ((v[0] * vold[0] + v[1] * vold[1]) < 0.0)
	{
		direction *= -1.0;
	}

	p[0] = p[0] + direction * this->IntegrationStep * v[0];
	p[1] = p[1] + direction * this->IntegrationStep * v[1];
	p[2] = p[2] + direction * this->IntegrationStep * v[2];
	points->InsertNextPoint(p[0], p[1], p[2]);

	vtkMath::Assign(v, vold);

	return direction;
}

void vtkExtractMidsurface3D::InsertNextCell(vtkCellArray *lines, vtkIdType id1, vtkIdType id2)
{
	vtkNew<vtkLine> line;
	line->GetPointIds()->SetId(0, id1);
	line->GetPointIds()->SetId(1, id2);
	lines->InsertNextCell(line);
}

int vtkExtractMidsurface3D::AppendPoints(double *p, int k, double direction, vtkImageData *input, vtkPoints *points, vtkCellArray *lines, vtkDataArray *tensArr, vtkDataArray *val)
{
	vtkIdType id = input->FindPoint(p);
	double *t = tensArr->GetTuple9(id);

	double n[3];  // eigenvector normal to surface
	double v1[3]; // eigenvector in direction of eigenvector with smallest eigenvalue
	double v2[3]; // eigenvector in direction of eigenvector with medium eigenvalue

	ComputeEigenvectors(v1, v2, n, tensArr, id);

	double vold[3];
	vtkMath::Assign(v2, vold);

	direction = InsertNextPoint(p, v2, direction, vold, points);
	InsertNextCell(lines, 0, k);

	double curval = val->GetTuple1(id);

	// stop when outside of segmentation mask
	while (curval > 0.5)
	{
		k++;

		double v[3];
		id = input->FindPoint(p);
		ComputeEigenvectors(v1, v2, n, tensArr, id);
		vtkMath::MultiplyScalar(n, -vtkMath::Dot(vold, n));
		vtkMath::Add(vold, n, v);

		direction = InsertNextPoint(p, v, direction, vold, points);
		InsertNextCell(lines, k - 1, k);

		curval = val->GetTuple1(id);
	}

	return k + 1;
}

void vtkExtractMidsurface3D::ComputeEigenvectors(double *v1, double *v2, double *n, vtkDataArray *tensArr, int k)
{
	double e[3]; // eigenvalues, not needed so far

	this->a[0][0] = tensArr->GetTuple9(k)[0];
	this->a[0][1] = tensArr->GetTuple9(k)[1];
	this->a[0][2] = tensArr->GetTuple9(k)[2];
	this->a[1][0] = tensArr->GetTuple9(k)[3];
	this->a[1][1] = tensArr->GetTuple9(k)[4];
	this->a[1][2] = tensArr->GetTuple9(k)[5];
	this->a[2][0] = tensArr->GetTuple9(k)[6];
	this->a[2][1] = tensArr->GetTuple9(k)[7];
	this->a[2][2] = tensArr->GetTuple9(k)[8];

	vtkMath::Jacobi(this->a, e, this->w);

	v1[0] = this->w[0][0];
	v1[1] = this->w[1][0];
	v1[2] = this->w[2][0];
	vtkMath::Normalize(v1);

	v2[0] = this->w[0][1];
	v2[1] = this->w[1][1];
	v2[2] = this->w[2][1];
	vtkMath::Normalize(v2);

	n[0] = this->w[0][2];
	n[1] = this->w[1][2];
	n[2] = this->w[2][2];
	vtkMath::Normalize(n);
}