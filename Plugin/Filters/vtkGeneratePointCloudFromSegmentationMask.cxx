#include "vtkGeneratePointCloudFromSegmentationMask.h"
#include <vtkImageData.h>
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

vtkStandardNewMacro(vtkGeneratePointCloudFromSegmentationMask);

//----------------------------------------------------------------------------
vtkGeneratePointCloudFromSegmentationMask::vtkGeneratePointCloudFromSegmentationMask() : InputArray(nullptr)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkGeneratePointCloudFromSegmentationMask::vtkGeneratePointCloudFromSegmentationMask() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), ":vtkGeneratePointCloudFromSegmentationMask:vtkGeneratePointCloudFromSegmentationMask() END");
}

vtkGeneratePointCloudFromSegmentationMask::~vtkGeneratePointCloudFromSegmentationMask()
{
	this->SetInputArray(nullptr);
}

int vtkGeneratePointCloudFromSegmentationMask::FillOutputPortInformation(int, vtkInformation *info)
{
	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
	return 1;
}

int vtkGeneratePointCloudFromSegmentationMask::RequestData(vtkInformation *vtkNotUsed(request),
														   vtkInformationVector **inputVector,
														   vtkInformationVector *outputVector)
{
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	auto input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	auto output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkNew<vtkPoints> points;
	vtkNew<vtkCellArray> cells;

	double p[3];
	vtkIdType num = 0;

	auto inArr = input->GetPointData()->GetArray(this->InputArray);

	for (auto i = 0; i < input->GetNumberOfPoints(); i++)
	{
		if (inArr->GetTuple1(i) > 0.5)
		{
			input->GetPoint(i, p);
			points->InsertNextPoint(p);
			vtkIdType pid[1];
			pid[0] = num;
			cells->InsertNextCell(1, pid);
			num++;
		}
	}

	vtkNew<vtkPolyData> poly;
	poly->SetPoints(points);
	poly->SetLines(cells);

	output->ShallowCopy(poly);

	return 1;
}
