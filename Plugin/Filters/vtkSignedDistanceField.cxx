/**
 * https://github.com/fau-vc/optimal-frames
 */

#include "vtkSignedDistanceField.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <queue>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkImageCast.h>
#include <vtkImageConvolve.h>
#include <vtkImageData.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkLogger.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkType.h>
#include <vtkImageEuclideanDistance.h>
#include <vtkImageMathematics.h>
#include <vtkExtractVOI.h>
#include <vtkImageShiftScale.h>

vtkStandardNewMacro(vtkSignedDistanceField);

vtkSignedDistanceField::vtkSignedDistanceField()
	: InputArray(nullptr), FinalResultName(nullptr) {
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
}

vtkSignedDistanceField::~vtkSignedDistanceField() {
	if (this->InputArray) {
		this->InputArray = nullptr;
	}
	if (this->SDFResultName) {
		this->SDFResultName = nullptr;
	}
	if (this->FinalResultName) {
		this->FinalResultName = nullptr;
	}
}

int vtkSignedDistanceField::RequestData(vtkInformation* vtkNotUsed(request),
	vtkInformationVector** inputVector,
	vtkInformationVector* outputVector) {
	// Get the info objects
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	vtkImageData* input =
		vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkImageData* output =
		vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkNew<vtkImageData> image;
	image->DeepCopy(input); // for initialization

	// Get the information on the domain
	dims = image->GetDimensions();
	this->width = dims[0];
	this->height = dims[1];
	this->depth = dims[2];

	// Get the input array
	vtkDataArray* inArr = input->GetPointData()->GetArray(this->InputArray);

	// array for signed distance map
	vtkNew<vtkDoubleArray> dist;
	dist->SetNumberOfComponents(1);
	dist->SetNumberOfTuples(inArr->GetNumberOfTuples());

	if (this->DistanceType == _2D)
	{
		vtkNew<vtkImageAppend> distanceAppend;
		distanceAppend->SetAppendAxis(2);

		for (int z = 0; z < this->depth; z++) {

			vtkLog(INFO, "depth: " << z);

			vtkNew<vtkExtractVOI> slice;
			slice->SetInputData(image);
			slice->SetVOI(0, this->width - 1, 0, this->height - 1, z, z);
			slice->Update();

			ComputeSignedDistanceField(dist, slice->GetOutput(), distanceAppend);
			distanceAppend->Update();
		}

		dist->ShallowCopy(distanceAppend->GetOutput()->GetPointData()->GetScalars());
	}
	else if (this->DistanceType == _3D)
	{
		ComputeSignedDistanceField(dist, image);
		dist->SetName(this->SDFResultName);
	}
	else {
		vtkLog(ERROR, "Invalid field type");
		return 0;
	}

	image->GetPointData()->AddArray(dist);

	if (this->Smoothing)
	{
		auto smooth_image = vtkSmartPointer<vtkDoubleArray>::New();
		smooth_image->SetNumberOfComponents(1);
		smooth_image->SetNumberOfTuples(inArr->GetNumberOfTuples());

		GaussianSmooth(smooth_image, dist);
		smooth_image->SetName(this->FinalResultName);

		image->GetPointData()->AddArray(smooth_image);
	} else {
		dist->SetName(this->FinalResultName);
	}

	output->ShallowCopy(image);

	vtkLog(INFO, "Finished");
	return 1;
}

void vtkSignedDistanceField::GaussianSmooth(vtkDoubleArray* smooth_image, vtkDoubleArray* dist)
{
	auto dist_image = vtkSmartPointer<vtkImageData>::New();
	dist_image->SetDimensions(this->dims);
	dist_image->GetPointData()->SetScalars(dist);

	vtkLog(INFO, "Smoothing the distance field");

	if (!this->ManualSmoothing) {
		vtkLog(INFO, "Calculating the distance field smoothing parameters automatically");
		ComputeStandardDeviations(dist);
		ComputeRadiusFactors();
	}

	vtkLog(INFO, "Standard deviations: " << this->StandardDeviations[0] << " " << this->StandardDeviations[1] << " " << this->StandardDeviations[2]);
	vtkLog(INFO, "Radius factors: " << this->RadiusFactors[0] << " " << this->RadiusFactors[1] << " " << this->RadiusFactors[2]);

	auto smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
	smooth->SetInputData(dist_image);
	smooth->SetStandardDeviations(this->StandardDeviations);
	smooth->SetRadiusFactors(this->RadiusFactors);
	smooth->Update();

	if (this->Normalize) {
		// normalize the smooth image (is this necessary for the consecutive steps?)
		auto shiftscale = vtkSmartPointer<vtkImageShiftScale>::New();
		shiftscale->SetInputData(smooth->GetOutput());
		double oldRange = smooth->GetOutput()->GetScalarRange()[1] - smooth->GetOutput()->GetScalarRange()[0];
		shiftscale->SetShift(-1.0 * smooth->GetOutput()->GetScalarRange()[0]);
		shiftscale->SetScale(1.0 / oldRange);
		shiftscale->Update();
		smooth_image->ShallowCopy(shiftscale->GetOutput()->GetPointData()->GetScalars());
	}
	else {
		smooth_image->ShallowCopy(smooth->GetOutput()->GetPointData()->GetScalars());
	}
}

void vtkSignedDistanceField::ComputeStandardDeviations(vtkDoubleArray* dist)
{
	// get maximum from the distance field
	double max = 0;
	vtkLog(INFO, "Computing standard deviations");

	for (int i = 0; i < dist->GetNumberOfTuples(); i++)
	{
		double value = dist->GetValue(i);
		if (value > max)
		{
			max = value;
		}
	}

	if (max < 1)
	{
		this->StandardDeviations[0] = 1;
		this->StandardDeviations[1] = 1;
		this->StandardDeviations[2] = 1;
		return;
	}
	else {
		this->StandardDeviations[0] = max / this->SigmaDivisor;
		this->StandardDeviations[1] = max / this->SigmaDivisor;
		this->StandardDeviations[2] = max / this->SigmaDivisor;
	}
}
void vtkSignedDistanceField::ComputeRadiusFactors()
{
	// this does not make sense, it was here because of the old implementation of the gaussian smoothing, now that we use vtk one, 
	// this calculation is done in the background (when setting the radius, they take that as the diameter as opposed to the radius and 
	// divide square of sigma by the value*2 -1)
	// if (this->TestingRadius <= 1)
	// {
	// 	this->TestingRadius = 1.1;
	// }
	// auto radius = (this->StandardDeviations[0] * 2) / (this->TestingRadius);
	vtkLog(INFO, "Radius: " << this->TestingRadius);
	this->RadiusFactors[0] = this->TestingRadius;
	this->RadiusFactors[1] = this->TestingRadius;
	this->RadiusFactors[2] = this->TestingRadius;
}

void vtkSignedDistanceField::ComputeSignedDistanceField(vtkDoubleArray* dist, vtkImageData* originalImage, vtkImageAppend* distanceAppend)
{
	vtkLog(INFO, "Computing signed distance field");
	auto root_inside = vtkSmartPointer<vtkImageMathematics>::New();
	auto root_outside = vtkSmartPointer<vtkImageMathematics>::New();

	vtkLog(INFO, "original image dimensions: " << originalImage->GetDimensions()[0] << " " << originalImage->GetDimensions()[1] << " " << originalImage->GetDimensions()[2]);

	if (this->FieldType == INSIDE || this->FieldType == BOTH)
	{
		vtkLog(INFO, "debug");
		auto dist_inside = vtkSmartPointer<vtkImageEuclideanDistance>::New();
		dist_inside->SetInputData(originalImage);
		dist_inside->SetDimensionality(this->DistanceType == _2D ? 2 : 3);
		dist_inside->SetAlgorithmToSaito();
		dist_inside->Update();
		vtkLog(INFO, "debug");

		root_inside->SetInput1Data(dist_inside->GetOutput());
		root_inside->SetOperationToSquareRoot();
		root_inside->Update();
		vtkLog(INFO, "debug");

	}

	if (this->FieldType == OUTSIDE || this->FieldType == BOTH)
	{
		auto invert = vtkSmartPointer<vtkImageMathematics>::New();
		invert->SetInput1Data(originalImage);
		invert->SetOperationToAddConstant();
		invert->SetConstantC(1);
		invert->Update();

		auto replace = vtkSmartPointer<vtkImageMathematics>::New();
		replace->SetInputData(invert->GetOutput());
		replace->SetOperationToReplaceCByK();
		replace->SetConstantK(0);
		replace->SetConstantC(2);
		replace->Update();

		auto dist_outside = vtkSmartPointer<vtkImageEuclideanDistance>::New();
		dist_outside->SetInputConnection(replace->GetOutputPort());
		dist_outside->SetDimensionality(this->DistanceType == _2D ? 2 : 3);
		dist_outside->SetAlgorithmToSaito();
		dist_outside->Update();

		root_outside->SetInput1Data(dist_outside->GetOutput());
		root_outside->SetOperationToSquareRoot();
		root_outside->Update();
	}

	if (this->FieldType == BOTH)
	{
		auto combine = vtkSmartPointer<vtkImageMathematics>::New();
		combine->SetInput1Data(root_inside->GetOutput());
		combine->SetInput2Data(root_outside->GetOutput());
		combine->SetOperationToSubtract();
		combine->Update();

		if (this->DistanceType == _2D)
		{
			distanceAppend->AddInputData(combine->GetOutput());
		}
		else if (this->DistanceType == _3D)
		{
			dist->ShallowCopy(combine->GetOutput()->GetPointData()->GetScalars());
		}
	}
	else if (this->FieldType == INSIDE)
	{
		if (this->DistanceType == _2D)
		{
		vtkLog(INFO, "debug");

			distanceAppend->AddInputData(root_inside->GetOutput());
			vtkLog(INFO, "distance append dimensions: " << distanceAppend->GetOutput()->GetDimensions()[0] << " " << distanceAppend->GetOutput()->GetDimensions()[1] << " " << distanceAppend->GetOutput()->GetDimensions()[2]);
		}
		else if (this->DistanceType == _3D)
		{
			dist->ShallowCopy(root_inside->GetOutput()->GetPointData()->GetScalars());
		}
	}
	else if (this->FieldType == OUTSIDE)
	{
		if (this->DistanceType == _2D)
		{
			distanceAppend->AddInputData(root_outside->GetOutput());
		}
		else if (this->DistanceType == _3D)
		{
			dist->ShallowCopy(root_outside->GetOutput()->GetPointData()->GetScalars());
		}
	}
}

int vtkSignedDistanceField::RequestInformation(
	vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector,
	vtkInformationVector* outputVector) {
	// get the info objects
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

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