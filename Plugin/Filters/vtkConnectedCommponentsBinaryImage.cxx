#include "vtkConnectedCommponentsBinaryImage.h"
#include <vtkPVLogger.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>

vtkStandardNewMacro(vtkConnectedCommponentsBinaryImage);

//----------------------------------------------------------------------------
vtkConnectedCommponentsBinaryImage::vtkConnectedCommponentsBinaryImage() : InputArray(nullptr)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkConnectedCommponentsBinaryImage::vtkConnectedCommponentsBinaryImage() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkConnectedCommponentsBinaryImage::vtkConnectedCommponentsBinaryImage() END");
}

vtkConnectedCommponentsBinaryImage::~vtkConnectedCommponentsBinaryImage()
{
	this->SetInputArray(nullptr);
}

int vtkConnectedCommponentsBinaryImage::RequestData(vtkInformation *vtkNotUsed(request),
													vtkInformationVector **inputVector,
													vtkInformationVector *outputVector)
{
	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	auto input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	auto output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkLog(INFO, "vtkConnectedCommponentsBinaryImage::RequestData() BEGIN");

	vtkLog(INFO, "InputArray: " << this->InputArray);

	output->DeepCopy(input);

	vtkLog(INFO, "Adding region IDs to the binary image.");

	AddRegionIDs(output); // requires the original input array
	vtkLog(INFO, "Region IDs added to the binary image.");
	return 1;
}

/// @brief from https://stackoverflow.com/questions/22051069/how-do-i-find-the-connected-components-in-a-binary-image
/// @param output the vtkImageData to add the RegionsIDs to
void vtkConnectedCommponentsBinaryImage::AddRegionIDs(vtkImageData *output)
{
	vtkLog(INFO, "1");
	int dims[3];
	output->GetDimensions(dims);

	vtkLog(INFO, "2");
	std::vector<std::vector<int>> g;
	g.resize(dims[0] + 2, std::vector<int>(dims[1] + 2, 0)); // add zero boundary

	vtkLog(INFO, "3");
	auto mask = output->GetPointData()->GetArray(this->InputArray);

	vtkLog(INFO, "4");
	for (int j = 0; j < dims[1]; j++)
		for (int i = 0; i < dims[0]; i++)
			g[i + 1][j + 1] = mask->GetTuple1(j * dims[0] + i);

	std::vector<std::vector<int>> w;
	w.resize(dims[0] + 2, std::vector<int>(dims[1] + 2, 0)); // add zero boundary
	
	vtkLog(INFO, "5");
	int set = 1;

	for (int i = 1; i <= dims[0]; i++)
		for (int j = 1; j <= dims[1]; j++)
			if (g[i][j] && !w[i][j])
				dfs(i, j, set++, g, w);

	vtkLog(INFO, "6");
	this->NumberOfConnectedRegions = set - 1;
	
	vtkNew<vtkIntArray> cc;
	cc->SetName("RegionID");
	cc->SetNumberOfComponents(1);
	cc->SetNumberOfTuples(mask->GetNumberOfTuples());

	vtkLog(INFO, "7");
	for (int j = 0; j < dims[1]; j++)
		for (int i = 0; i < dims[0]; i++)
			cc->SetTuple1(j * dims[0] + i, w[i + 1][j + 1]);

	vtkLog(INFO, "8");
	output->GetPointData()->AddArray(cc);
	output->GetAttributes(vtkDataObject::POINT)->SetActiveScalars("RegionID");
}

void vtkConnectedCommponentsBinaryImage::dfs(int x, int y, int c, std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &w)
{
	// first four for 4-connectivity, all for 8-connectivity
	int dx[8] = {0, 1, 0, -1, 1, 1, -1, -1};
	int dy[8] = {1, 0, -1, 0, 1, -1, -1, 1};

	w[x][y] = c;
	for (unsigned int i = 0; i < this->Connectivity; i++)
	{
		int nx = x + dx[i], ny = y + dy[i];
		if ((g[nx][ny] > 0) && (w[nx][ny] == 0))
			dfs(nx, ny, c, g, w);
	}
}

void vtkConnectedCommponentsBinaryImage::DilateConnectedComponents(vtkImageData* slice, vtkIntArray* dilatedMaskArr)
{
	auto mask = slice->GetPointData()->GetArray("RegionID");
	vtkLog(INFO, "Dilating the mask.");
	auto size = mask->GetNumberOfTuples();
	auto width = slice->GetDimensions()[0];
	auto height = slice->GetDimensions()[1];
	for (int i = 0; i < size; i++)
	{	
		bool found = false;
		int reg = 0;
		auto pixel = mask->GetTuple1(i);
		if (pixel == 0)
		{
			int x = i % width;
			int y = i / width;
			for (int j = -1; j <= 1; j++)
			{
				for (int k = -1; k <= 1; k++)
				{
					int nx = x + k;
					int ny = y + j;
					if (nx >= 0 && nx < width && ny >= 0 && ny < height) 
					{
						int idx = ny * width + nx;
						auto val = mask->GetTuple1(idx);
						if (val > 0)
						{
							found = true;
							reg = val;
						} 
					}
				}
			}
			if (found) 
			{
				dilatedMaskArr->SetTuple1(i, reg);
			} else {
				dilatedMaskArr->SetTuple1(i, 0);
			}
		} else 
		{
			dilatedMaskArr->SetTuple1(i, pixel);
		}
	}
}

void vtkConnectedCommponentsBinaryImage::CloseConnectedComponents(vtkImageData* slice, vtkIntArray* closedMaskArr)
{
	auto mask = slice->GetPointData()->GetArray("dilatedMask");
	vtkLog(INFO, "Closing the mask.");
	auto size = mask->GetNumberOfTuples();
	auto width = slice->GetDimensions()[0];
	auto height = slice->GetDimensions()[1];
	for (int i = 0; i < size; i++)
	{	
		auto pixel = mask->GetTuple1(i);

		if (pixel > 0) { // Only process non-zero pixels for erosion
        int x = i % width;
        int y = i / width;
        bool edgePixel = false;
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                int nx = x + k;
                int ny = y + j;
                if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                    int idx = ny * width + nx;
                    auto val = mask->GetTuple1(idx);
                    if (val == 0) {
                        edgePixel = true;
                        break;
                    }
                }
            }
            if (edgePixel) {
                break;
            }
        }
        if (edgePixel) {
            // Erode pixel by setting its value to zero
            closedMaskArr->SetTuple1(i, 0);
        } else 
		{
			closedMaskArr->SetTuple1(i, pixel);
		}
    }
	else 
	{
		closedMaskArr->SetTuple1(i, pixel);
	}
	}
}