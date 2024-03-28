#ifndef __vtkZipperTriangulation_h
#define __vtkZipperTriangulation_h
 
#include "vtkMidsurfaceExtractorModule.h" // for export
#include <vtkPVLogger.h>

#include "vtkPolyDataAlgorithm.h"
 
class VTKMIDSURFACEEXTRACTOR_EXPORT vtkZipperTriangulation : public vtkPolyDataAlgorithm
//class STTKPLUGIN_EXPORT vtkExtractCoreLines : public vtkImageAlgorithm
{
public:
	static vtkZipperTriangulation *New();
	vtkTypeMacro(vtkZipperTriangulation, vtkPolyDataAlgorithm);

	vtkSetMacro(ZipperAlpha, double);
	vtkGetMacro(ZipperAlpha, double);

protected:

	vtkZipperTriangulation();
	~vtkZipperTriangulation();

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
	int FillInputPortInformation(int port, vtkInformation *info) override;

private:
	vtkZipperTriangulation(const vtkZipperTriangulation&) = delete;
	void operator=(const vtkZipperTriangulation&) = delete;

	double ZipperAlpha;
};
 
#endif // __vtkZipperTriangulation_h
