#include "vtkZipperTriangulation.h"
#include <vtkPolyData.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

vtkStandardNewMacro(vtkZipperTriangulation);

vtkZipperTriangulation::vtkZipperTriangulation()
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkZipperTriangulation::vtkZipperTriangulation() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), ":vtkZipperTriangulation:vtkZipperTriangulation() END");
}

vtkZipperTriangulation::~vtkZipperTriangulation()
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkZipperTriangulation::~vtkZipperTriangulation() BEGIN");

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), ":vtkZipperTriangulation:~vtkZipperTriangulation() END");
}

int vtkZipperTriangulation::RequestData(vtkInformation *vtkNotUsed(request),
										vtkInformationVector **inputVector,
										vtkInformationVector *outputVector)
{
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	auto input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	auto output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	output->ShallowCopy(input);

	return 1;
}

int vtkZipperTriangulation::FillInputPortInformation(int port, vtkInformation *info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}
