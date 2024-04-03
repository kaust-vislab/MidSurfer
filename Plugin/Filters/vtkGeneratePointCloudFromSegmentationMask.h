#ifndef __vtkGeneratePointCloudFromSegmentationMask_h
#define __vtkGeneratePointCloudFromSegmentationMask_h

#include "vtkMidSurferModule.h" // for export
#include <vtkPVLogger.h>

#include "vtkImageAlgorithm.h"
#include <vtkPoints.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkAppendPolyData.h>
#include <vtkConnectivityFilter.h>

class VTKMIDSURFER_EXPORT vtkGeneratePointCloudFromSegmentationMask : public vtkImageAlgorithm
{
public:
	static vtkGeneratePointCloudFromSegmentationMask *New();
	vtkTypeMacro(vtkGeneratePointCloudFromSegmentationMask, vtkImageAlgorithm);

	vtkGetStringMacro(InputArray);
	vtkSetStringMacro(InputArray);
	
	vtkSetMacro(ComputeNormals, bool);
	vtkGetMacro(ComputeNormals, bool);

protected:
	vtkGeneratePointCloudFromSegmentationMask();
	~vtkGeneratePointCloudFromSegmentationMask();

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
	int FillOutputPortInformation(int, vtkInformation *) override;

private:
	vtkGeneratePointCloudFromSegmentationMask(const vtkGeneratePointCloudFromSegmentationMask &) = delete;
	void operator=(const vtkGeneratePointCloudFromSegmentationMask &) = delete;

	// GUI parameters
	char *InputArray;
	bool ComputeNormals;
};

#endif // __vtkGeneratePointCloudFromSegmentationMask_h
