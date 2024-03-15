#include "vtkConnectedCommponentsBinaryImage3D.h"
#include <vtkPVLogger.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>

vtkStandardNewMacro(vtkConnectedCommponentsBinaryImage3D);

//----------------------------------------------------------------------------
vtkConnectedCommponentsBinaryImage3D::vtkConnectedCommponentsBinaryImage3D() : InputArray(nullptr)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkConnectedCommponentsBinaryImage3D::vtkConnectedCommponentsBinaryImage3D() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkConnectedCommponentsBinaryImage3D::vtkConnectedCommponentsBinaryImage3D() END");
}

vtkConnectedCommponentsBinaryImage3D::~vtkConnectedCommponentsBinaryImage3D()
{
	this->SetInputArray(nullptr);
}

int vtkConnectedCommponentsBinaryImage3D::RequestData(vtkInformation *vtkNotUsed(request),
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

	AddRegionIDs(output); // requires the original input array

	return 1;
}

/// @brief from https://stackoverflow.com/questions/22051069/how-do-i-find-the-connected-components-in-a-binary-image
/// @param output the vtkImageData to add the RegionsIDs to
void vtkConnectedCommponentsBinaryImage3D::AddRegionIDs(vtkImageData *output)
{
	int dims[3];
	output->GetDimensions(dims);

	std::vector<std::vector<std::vector<int>>> g;
	g.resize(dims[0] + 2, std::vector<std::vector<int>>(dims[1] + 2, std::vector<int>(dims[2] + 2, 0)));

	auto mask = output->GetPointData()->GetArray(this->InputArray);

	for (int k = 0; k < dims[2]; k++)
		for (int j = 0; j < dims[1]; j++)
			for (int i = 0; i < dims[0]; i++)
				g[i + 1][j + 1][k + 1] = mask->GetTuple1(k * dims[0] * dims[1] + j * dims[0] + i);

	std::vector<std::vector<std::vector<int>>> w;
	w.resize(dims[0] + 2, std::vector<std::vector<int>>(dims[1] + 2, std::vector<int>(dims[2] + 2, 0)));

	int set = 1;

	for (int i = 1; i <= dims[0]; i++)
		for (int j = 1; j <= dims[1]; j++)
			for (int k = 1; k <= dims[2]; k++)
				if ((g[i][j][k] > 0) && (w[i][j][k] == 0))
					dfs(i, j, k, set++, g, w);

	this->NumberOfConnectedRegions = set - 1;

	vtkNew<vtkIntArray> cc;
	cc->SetName("RegionID");
	cc->SetNumberOfComponents(1);
	cc->SetNumberOfTuples(mask->GetNumberOfTuples());

	for (int k = 0; k < dims[2]; k++)
		for (int j = 0; j < dims[1]; j++)
			for (int i = 0; i < dims[0]; i++)
				cc->SetTuple1(k * dims[0] * dims[1] + j * dims[0] + i, w[i + 1][j + 1][k + 1]);

	output->GetPointData()->AddArray(cc);
	output->GetAttributes(vtkDataObject::POINT)->SetActiveScalars("RegionID");
}

void vtkConnectedCommponentsBinaryImage3D::dfs(int x, int y, int z, int c, std::vector<std::vector<std::vector<int>>> &g, std::vector<std::vector<std::vector<int>>> &w)
{
	// first four for 4-connectivity, all for 8-connectivity
	int dx[26] = {1, 0, 0, -1, 0, 0, /**/ 0, 0, 0, 0, /**/ 1, 1, -1, -1, /**/ 1, 1, -1, -1, /**/ 1, 1, 1, 1, -1, -1, -1, -1};
	int dy[26] = {0, 1, 0, 0, -1, 0, /**/ 1, 1, -1, -1, /**/ 0, 0, 0, 0, /**/ 1, -1, 1, -1, /**/ 1, 1, -1, -1, 1, 1, -1, -1};
	int dz[26] = {0, 0, 1, 0, 0, -1, /**/ 1, -1, 1, -1, /**/ 1, -1, 1, -1, /**/ 0, 0, 0, 0, /**/ 1, -1, 1, -1, 1, -1, 1, -1};

	w[x][y][z] = c;
	for (int i = 0; i < this->Connectivity; i++)
	{
		int nx = x + dx[i], ny = y + dy[i], nz = z + dz[i];
		if ((g[nx][ny][nz] > 0) && (w[nx][ny][nz] == 0))
			dfs(nx, ny, nz, c, g, w);
	}
}
