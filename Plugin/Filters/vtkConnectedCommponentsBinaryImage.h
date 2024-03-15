#ifndef __vtkConnectedCommponentsBinaryImage_h
#define __vtkConnectedCommponentsBinaryImage_h

#include "vtkMidsurfaceExtractorModule.h" // for export
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include "vtkImageAlgorithm.h"

class VTKMIDSURFACEEXTRACTOR_EXPORT vtkConnectedCommponentsBinaryImage : public vtkImageAlgorithm
{
public:
	static vtkConnectedCommponentsBinaryImage *New();
	vtkTypeMacro(vtkConnectedCommponentsBinaryImage, vtkImageAlgorithm);

	vtkGetStringMacro(InputArray);
	vtkSetStringMacro(InputArray);

	vtkSetMacro(Connectivity, unsigned int);
	vtkGetMacro(Connectivity, unsigned int);

	unsigned int GetNumberOfConnectedRegions() { return this->NumberOfConnectedRegions; }
	void DilateConnectedComponents(vtkImageData* slice, vtkIntArray* dilatedMaskArr);
	void CloseConnectedComponents(vtkImageData* slice, vtkIntArray* closedMaskArr);

protected:
	vtkConnectedCommponentsBinaryImage();
	~vtkConnectedCommponentsBinaryImage();

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
	vtkConnectedCommponentsBinaryImage(const vtkConnectedCommponentsBinaryImage &) = delete;
	void operator=(const vtkConnectedCommponentsBinaryImage &) = delete;

	void AddRegionIDs(vtkImageData *output);
	void dfs(int x, int y, int c, std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &w);

	unsigned int NumberOfConnectedRegions;

	// GUI paramters
	char *InputArray;
	unsigned int Connectivity;
};

#endif // __vtkConnectedCommponentsBinaryImage_h
