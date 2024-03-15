#ifndef __vtkExtractMidsurface3D_h
#define __vtkExtractMidsurface3D_h

#include "vtkMidsurfaceExtractorModule.h" // for export
#include <vtkPVLogger.h>

#include "vtkImageAlgorithm.h"
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkAppendPolyData.h>

class VTKMIDSURFACEEXTRACTOR_EXPORT vtkExtractMidsurface3D : public vtkImageAlgorithm
{
public:
	static vtkExtractMidsurface3D *New();
	vtkTypeMacro(vtkExtractMidsurface3D, vtkImageAlgorithm);

	vtkGetStringMacro(InputArray);
	vtkSetStringMacro(InputArray);

	vtkSetMacro(IntegrationStep, double);
	vtkGetMacro(IntegrationStep, double);

	vtkSetVector3Macro(RadiusFactors, double);
	vtkGetVector3Macro(RadiusFactors, double);

	vtkSetVector3Macro(StandardDeviations, double);
	vtkGetVector3Macro(StandardDeviations, double);

protected:
	vtkExtractMidsurface3D();
	~vtkExtractMidsurface3D();

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
	int FillOutputPortInformation(int, vtkInformation *) override;

private:
	vtkExtractMidsurface3D(const vtkExtractMidsurface3D &) = delete;
	void operator=(const vtkExtractMidsurface3D &) = delete;

	template <class TReal>
	TReal **create_matrix(long nrow, long ncol);
	template <class TReal>
	void free_matrix(TReal **m);
	double InsertNextPoint(double *p, double *v, double direction, double *vold, vtkPoints *points);
	void InsertNextCell(vtkCellArray *lines, vtkIdType id1, vtkIdType id2);
	int AppendPoints(double *p, int k, double direction, vtkImageData *input, vtkPoints *points, vtkCellArray *lines, vtkDataArray *tensArr, vtkDataArray *val);
	vtkIdType FindMaximumLocation(vtkDataArray *arr);
	void ExtractMidsurface(vtkImageData *input, vtkAppendPolyData *append);
	void ComputeEigenvectors(double *v1, double *v2, double *n, vtkDataArray *tensArr, int k);

	std::vector<std::array<double, 3>> StartPoints;
	double **a; // matrix for tensor
	double **w; // matrix for eigenvectors of a 

	char *InputArray;
	double IntegrationStep;
	double RadiusFactors[3];
	double StandardDeviations[3];
};

#endif // __vtkExtractMidsurface3D_h
