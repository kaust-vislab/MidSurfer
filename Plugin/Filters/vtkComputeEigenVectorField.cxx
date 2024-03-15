#include "vtkComputeEigenVectorField.h"
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkVector.h>
#include <vtkMath.h>
#include <vtkImageCast.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkContourFilter.h>
#include <vtkPointLocator.h>
#include <vtkPassSelectedArrays.h>
#include <vtkDataArraySelection.h>

#include "vtkSignedDistanceField.h"

vtkStandardNewMacro(vtkComputeEigenVectorField);

//----------------------------------------------------------------------------
vtkComputeEigenVectorField::vtkComputeEigenVectorField() : InputArray(nullptr)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkComputeEigenVectorField::vtkComputeEigenVectorField() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkComputeEigenVectorField::vtkComputeEigenVectorField() END");
}

vtkComputeEigenVectorField::~vtkComputeEigenVectorField()
{
	this->SetInputArray(nullptr);
}

template <class TReal>
TReal **vtkComputeEigenVectorField::create_matrix(long nrow, long ncol)
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
void vtkComputeEigenVectorField::free_matrix(TReal **m)
{
	free(m[0]);
	delete[] m;
}

void  vtkComputeEigenVectorField::dfs(int x, int y, int c, std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &w)
{
	int dx[8] = {-1, 0, 1, 1, 1, 0, -1, -1};
	int dy[8] = {1, 1, 1, 0, -1, -1, -1, 0};

	w[x][y] = c;
	for (int i = 0; i < 8; i++)
	{
		int nx = x + dx[i], ny = y + dy[i];
		if (g[nx][ny] && !w[nx][ny])
			dfs(nx, ny, c, g, w);
	}
}

int vtkComputeEigenVectorField::RequestData(vtkInformation *vtkNotUsed(request),
											vtkInformationVector **inputVector,
											vtkInformationVector *outputVector)
{
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	auto input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	auto output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	output->DeepCopy(input);

    // Compute gradients
	auto grad = vtkSmartPointer<vtkGradientFilter>::New();
	if (this->SmoothInput == ESmoothInput::SMOOTH_INPUT_OFF)
	{
		grad->SetInputData(output);
	}
	else if (this->SmoothInput == ESmoothInput::SMOOTH_INPUT_GAUSSIAN)
	{
		ComputeGaussianSmoothing(output, grad);
	}
	else if (this->SmoothInput == ESmoothInput::SMOOTH_INPUT_SIGNED_DISTANCE_FIELD_SIMPLE)
	{
		ComputeSimpleSignedDistanceField(output, grad);
	}
	else if (this->SmoothInput == ESmoothInput::SMOOTH_INPUT_SIGNED_DISTANCE_FIELD_MINKOWSI)
	{
		ComputeSignedDistanceField(output, grad);
	}
	else
	{
		vtkErrorMacro("Unknown smoothing method");
		return 0;
	}

	grad->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, this->InputArray);
	grad->SetResultArrayName("Gradient");
	grad->Update();

	// Compute tensors
	auto tens = vtkSmartPointer<vtkGradientFilter>::New();
	tens->SetInputConnection(grad->GetOutputPort());
	tens->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, "Gradient");
	tens->SetResultArrayName("T");
	tens->Update();

	auto pointData = tens->GetOutput()->GetPointData();
	auto gradArr = pointData->GetArray("Gradient");
	auto tensArr = pointData->GetArray("T");

	ComputeVectorField(gradArr, tensArr);

    AddRegionIDs(output); // requires the original input array
	output->GetPointData()->AddArray(pointData->GetArray(this->InputArray)); // replaces input array with smoothed one
	output->GetPointData()->AddArray(gradArr);
	output->GetPointData()->AddArray(tensArr);

	return 1;
}

/// @brief from https://stackoverflow.com/questions/22051069/how-do-i-find-the-connected-components-in-a-binary-image
/// @param output the vtkImageData to add the RegionsIDs to
void vtkComputeEigenVectorField::AddRegionIDs(vtkImageData *output)
{
    int dims[3];
    output->GetDimensions(dims);

    std::vector<std::vector<int>> g;
    g.resize(dims[0], std::vector<int>(dims[1], 0));

    auto mask = output->GetPointData()->GetArray(this->InputArray);

    for (int j = 0; j < dims[1]; j++)
        for (int i = 0; i < dims[0]; i++)
            g[i][j] = mask->GetTuple1(j * dims[0] + i);

    std::vector<std::vector<int>> w;
    w.resize(dims[0], std::vector<int>(dims[1], 0));

    int set = 1;

    for (int i = 1; i < dims[0] - 1; i++)
        for (int j = 1; j < dims[1] - 1; j++)
            if (g[i][j] && !w[i][j])
                dfs(i, j, set++, g, w);

    vtkNew<vtkIntArray> cc;
    cc->SetName("RegionID");
    cc->SetNumberOfComponents(1);
    cc->SetNumberOfTuples(mask->GetNumberOfTuples());

    for (int j = 0; j < dims[1]; j++)
        for (int i = 0; i < dims[0]; i++)
            cc->SetTuple1(j * dims[0] + i, w[i][j]);

    output->GetPointData()->AddArray(cc);
}

void vtkComputeEigenVectorField::ComputeSignedDistanceField(vtkImageData *output, vtkGradientFilter *grad)
{
	vtkNew<vtkSignedDistanceField> dist;
	dist->SetInputData(output);
	dist->SetInputArray(this->InputArray);
	dist->SetDistanceType(this->DistanceType);
	dist->SetFieldType(this->FieldType);
	dist->SetSmoothing(1);
	dist->SetNormalize(1);
	dist->SetSDFResultName("dist");
	dist->SetFinalResultName("smooth");
	dist->Update();

	vtkNew<vtkPassSelectedArrays> pass;
	pass->SetInputConnection(dist->GetOutputPort());
	pass->GetPointDataArraySelection()->EnableArray("dist");
	pass->Update();

	vtkImageData::SafeDownCast(pass->GetOutput())->GetPointData()->GetArray("dist")->SetName(this->InputArray);
	grad->SetInputData(pass->GetOutput());
}

void vtkComputeEigenVectorField::ComputeGaussianSmoothing(vtkImageData *output, vtkGradientFilter *grad)
{
	vtkNew<vtkImageCast> castFilter;
	castFilter->SetOutputScalarTypeToDouble();
	castFilter->SetInputData(output);
	castFilter->Update();

	vtkNew<vtkImageGaussianSmooth> smooth;
	smooth->SetInputData(castFilter->GetOutput());
	smooth->SetRadiusFactors(this->RadiusFactors);
	smooth->SetStandardDeviations(this->StandardDeviations);
	smooth->Update();

	grad->SetInputData(smooth->GetOutput());
}

void vtkComputeEigenVectorField::ComputeVectorField(vtkDataArray *gradArr, vtkDataArray *tensArr)
{
	double **a = create_matrix<double>(3, 3);
	double **w = create_matrix<double>(3, 3);

	for (int k = 0; k < gradArr->GetNumberOfTuples(); k++)
	{
		auto T = tensArr->GetTuple9(k);
		double v[3];
		a[0][0] = T[0];
		a[0][1] = T[1];
		a[0][2] = T[2];
		a[1][0] = T[3];
		a[1][1] = T[4];
		a[1][2] = T[5];
		a[2][0] = T[6];
		a[2][1] = T[7];
		a[2][2] = T[8];

		if (vtkMath::Jacobi(a, v, w))
		{
			if (w[2][0] > 0.0)
				gradArr->SetTuple3(k, w[0][1], w[1][1], w[2][1]);
			else
				gradArr->SetTuple3(k, w[0][0], w[1][0], w[2][0]);

			gradArr->SetName("vec");
		}
		else
		{
			vtkLog(INFO, "vtkMath::Jacobi returned false, a: ");
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					vtkLog(INFO, "a[" << std::to_string(i) << "][" << std::to_string(j) << "] = " << std::to_string(a[i][j]));
			gradArr->SetTuple3(k, 1.0, 0.0, 0.0); // put some default vector, should not happen anyway
		}
	}

	free_matrix(w);
	free_matrix(a);
}

void vtkComputeEigenVectorField::ComputeSimpleSignedDistanceField(vtkImageData *output, vtkGradientFilter *grad)
{
	vtkNew<vtkImageCast> castFilter;
	castFilter->SetOutputScalarTypeToDouble();
	castFilter->SetInputData(output);
	castFilter->Update();

	vtkNew<vtkImageGaussianSmooth> smooth;
	smooth->SetInputData(castFilter->GetOutput());
	smooth->SetRadiusFactors(this->RadiusFactorsSDF);
	smooth->SetStandardDeviations(this->StandardDeviationsSDF);
	smooth->Update();

	// Create the isosurface
	vtkNew<vtkContourFilter> contour;
	contour->SetInputConnection(smooth->GetOutputPort());
	contour->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, this->InputArray);
	contour->GenerateValues(1, this->Threshold, this->Threshold);
	contour->ComputeNormalsOff();
	contour->ComputeGradientsOff();
	contour->Update();

	// Initialize the locator
	vtkNew<vtkPointLocator> pointTree;
	pointTree->SetDataSet(contour->GetOutput());
	pointTree->BuildLocator();

	auto inArr = smooth->GetOutput()->GetPointData()->GetArray(this->InputArray);
	std::array<double, 3> p;  // for voxel location
	std::array<double, 3> pc; // for clostest point of contour

	double min = VTK_DOUBLE_MAX;
	double max = VTK_DOUBLE_MIN;

	for (auto i = 0; i < inArr->GetNumberOfTuples(); i++)
	{
		output->GetPoint(i, p.data());
		auto val = inArr->GetTuple1(i);

		contour->GetOutput()->GetPoint(pointTree->FindClosestPoint(p.data()), pc.data());
		double udist = std::sqrt(vtkMath::Distance2BetweenPoints(p.data(), pc.data()));
		udist = val < this->Threshold ? -1.0 * udist : udist;
		inArr->SetTuple1(i, udist);

		if (udist < min)
			min = udist;

		if (udist > max)
			max = udist;
	}

	for (auto i = 0; i < inArr->GetNumberOfTuples(); i++)
	{
		auto val = inArr->GetTuple1(i);

		if (val >= 0.0)
			inArr->SetTuple1(i, max * std::pow(val / max, this->ExponentDist));
		else
			inArr->SetTuple1(i, min * std::pow(val / min, this->ExponentDist));
	}

	grad->SetInputData(smooth->GetOutput());
}
