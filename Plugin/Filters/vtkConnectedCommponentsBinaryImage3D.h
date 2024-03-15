#ifndef __vtkConnectedCommponentsBinaryImage3D_h
#define __vtkConnectedCommponentsBinaryImage3D_h

#include "vtkMidsurfaceExtractorModule.h" // for export
#include <vtkImageData.h>

#include "vtkImageAlgorithm.h"

class VTKMIDSURFACEEXTRACTOR_EXPORT vtkConnectedCommponentsBinaryImage3D : public vtkImageAlgorithm
{
public:
	static vtkConnectedCommponentsBinaryImage3D *New();
	vtkTypeMacro(vtkConnectedCommponentsBinaryImage3D, vtkImageAlgorithm);

	vtkGetStringMacro(InputArray);
	vtkSetStringMacro(InputArray);

	vtkSetMacro(Connectivity, unsigned int);
	vtkGetMacro(Connectivity, unsigned int);

	unsigned int GetNumberOfCennectedRegions() { return this->NumberOfConnectedRegions; }

protected:
	vtkConnectedCommponentsBinaryImage3D();
	~vtkConnectedCommponentsBinaryImage3D();

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
	vtkConnectedCommponentsBinaryImage3D(const vtkConnectedCommponentsBinaryImage3D &) = delete;
	void operator=(const vtkConnectedCommponentsBinaryImage3D &) = delete;

	void AddRegionIDs(vtkImageData *output);
	void dfs(int x, int y, int nz, int c, std::vector<std::vector<std::vector<int>>> &g, std::vector<std::vector<std::vector<int>>> &w);

	unsigned int NumberOfConnectedRegions;

	// GUI paramters
	char *InputArray;
	unsigned int Connectivity;
};

#endif // __vtkConnectedCommponentsBinaryImage3D_h
