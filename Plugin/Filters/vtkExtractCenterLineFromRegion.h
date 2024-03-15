#ifndef __vtkExtractCenterLineFromRegion_h
#define __vtkExtractCenterLineFromRegion_h

#include "vtkMidsurfaceExtractorModule.h" // for export
#include <vtkPVLogger.h>

#include "vtkImageAlgorithm.h"
#include <vtkPoints.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkAppendPolyData.h>
#include <vtkConnectivityFilter.h>
#include <vtkImageInterpolator.h>
#include <vtkSmartPointer.h>

typedef std::array<double, 3> Point3;
typedef std::array<double, 3> Vector3;

class VTKMIDSURFACEEXTRACTOR_EXPORT vtkExtractCenterLineFromRegion : public vtkImageAlgorithm
{
public:
	enum ESmoothingMethod
	{
		SMOOTH_INPUT_OFF = 1,
		SMOOTH_INPUT_GAUSSIAN_2D = 2,
		SMOOTH_INPUT_GAUSSIAN_3D = 3,
		SMOOTH_INPUT_SIGNED_DISTANCE_FIELD_SIMPLE = 4,
		SMOOTH_INPUT_SIGNED_DISTANCE_FIELD = 5
	};
	enum EMorphological { DILATION = 1, CLOSING = 2, NONE = 3 };

	static vtkExtractCenterLineFromRegion *New();
	vtkTypeMacro(vtkExtractCenterLineFromRegion, vtkImageAlgorithm);

	vtkGetStringMacro(InputArray);
	vtkSetStringMacro(InputArray);

	vtkSetMacro(Exponent, double);
	vtkGetMacro(Exponent, double);

	vtkSetMacro(Threshold, double);
	vtkGetMacro(Threshold, double);

	vtkSetMacro(IntegrationStep, double);
	vtkGetMacro(IntegrationStep, double);

	vtkSetMacro(AutomaticStepSize, bool);
	vtkGetMacro(AutomaticStepSize, bool);

	vtkSetVector3Macro(StartPoint, double);
	vtkGetVector3Macro(StartPoint, double);

	vtkSetMacro(GoldenSectionSearch, bool);
	vtkGetMacro(GoldenSectionSearch, bool);

	vtkSetMacro(Tolerance, double);
	vtkGetMacro(Tolerance, double);

	vtkSetMacro(ResultType, int);
	vtkGetMacro(ResultType, int);

	vtkSetMacro(Morphological, unsigned int);
	vtkGetMacro(Morphological, unsigned int);

protected:
	vtkExtractCenterLineFromRegion();
	~vtkExtractCenterLineFromRegion();

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
	int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;
	int FillOutputPortInformation(int, vtkInformation *) override;

private:
	vtkExtractCenterLineFromRegion(const vtkExtractCenterLineFromRegion &) = delete;
	void operator=(const vtkExtractCenterLineFromRegion &) = delete;

	template <class TReal>
	TReal **create_matrix(const long nrow, const long ncol);
	template <class TReal>
	void free_matrix(TReal **m);
	void ComputeNextPoint(Point3 &p, Vector3 &v, Vector3 &vold);
	void InsertNextCell(vtkCellArray *lines, const vtkIdType id1, const vtkIdType id2);
	int AppendPoints(Point3 &p, int k, vtkImageData *input, vtkPoints *points, vtkCellArray *lines);
	// void ExtractCenterlineFromSlice(vtkImageData *input, vtkAppendPolyData *append);
	void ExtractCenterlineFromRegion(const Point3 &arr, vtkImageData *input, vtkPolyData *centerline);
	void ComputeEigenvector(Vector3 &v, const vtkIdType id);
	void DoGoldenSectionSearch(Point3 &p, Vector3 &v);
	void FindBounds(const Point3 &p, const Vector3 &v, double *ax, double *cx, double (vtkExtractCenterLineFromRegion::*func)(const Point3 &));
	double golden(Point3 &p, Vector3 &v, double ax, double bx, double cx, double (vtkExtractCenterLineFromRegion::*func)(const Point3 &, const Vector3 &, double));
	double GetSmoothedValue(const Point3 &p);
	double GetSmoothedValueAlongLine(const Point3 &p, const Vector3 &v, double x);

	vtkDataArray *heightArr;
	vtkDataArray *smoothHeightArr;
	// vtkDataArray *gradArr;
	vtkDataArray *tensArr;
	std::vector<Point3> StartPoints;
	double **a;
	double **w;
	double CurrentDirection;
	vtkSmartPointer<vtkImageInterpolator> interp;

	// GUI parameters
	char *InputArray;
	double Exponent;
	double Threshold;
	double IntegrationStep;
	bool AutomaticStepSize;
	int ResultType;
	bool GoldenSectionSearch;
	double Tolerance;
	double StartPoint[3];
	unsigned int Morphological;
};

#endif // __vtkExtractCenterLineFromRegion_h
