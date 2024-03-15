#ifndef __vtkExtractMidsurface_h
#define __vtkExtractMidsurface_h

#include "vtkMidsurfaceExtractorModule.h" // for export
#include <vtkPVLogger.h>

#include "vtkImageAlgorithm.h"
#include <vtkPoints.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkAppendPolyData.h>
#include <vtkConnectivityFilter.h>
#include <vtkSmartPointer.h>

#include "vtkExtractCenterLine.h"

typedef std::array<double, 3> Point3;
typedef std::array<double, 3> Vector3;

class VTKMIDSURFACEEXTRACTOR_EXPORT vtkExtractMidsurface : public vtkImageAlgorithm
{
public:
	enum ESmoothingMethod
	{
		SMOOTH_INPUT_OFF = 1,
		SMOOTH_INPUT_GAUSSIAN_2D = 2,
		SMOOTH_INPUT_GAUSSIAN_3D = 3,
		SMOOTH_INPUT_SIGNED_DISTANCE_FIELD_SIMPLE = 4,
		SMOOTH_INPUT_SIGNED_DISTANCE_FIELD = 5,
		SMOOTH_INPUT_SDF_3D = 6
	};

	static vtkExtractMidsurface* New();
	vtkTypeMacro(vtkExtractMidsurface, vtkImageAlgorithm);

	vtkGetStringMacro(InputArray);
	vtkSetStringMacro(InputArray);

	vtkSetMacro(SmoothInput, unsigned int);
	vtkGetMacro(SmoothInput, unsigned int);

	vtkSetMacro(DistanceType, unsigned int);
	vtkGetMacro(DistanceType, unsigned int);

	vtkSetMacro(FieldType, unsigned int);
	vtkGetMacro(FieldType, unsigned int);

	vtkSetMacro(Threshold, double);
	vtkGetMacro(Threshold, double);

	vtkSetMacro(IntegrationStep, double);
	vtkGetMacro(IntegrationStep, double);

	vtkSetMacro(AutomaticStepSize, bool);
	vtkGetMacro(AutomaticStepSize, bool);

	vtkSetVector3Macro(RadiusFactors, double);
	vtkGetVector3Macro(RadiusFactors, double);

	vtkSetVector3Macro(StandardDeviations, double);
	vtkGetVector3Macro(StandardDeviations, double);

	vtkSetVector3Macro(RadiusFactorsSDF, double);
	vtkGetVector3Macro(RadiusFactorsSDF, double);

	vtkSetVector3Macro(StandardDeviationsSDF, double);
	vtkGetVector3Macro(StandardDeviationsSDF, double);

	vtkSetMacro(GoldenSectionSearch, bool);
	vtkGetMacro(GoldenSectionSearch, bool);

	vtkSetMacro(ShapeDetection, bool);
	vtkGetMacro(ShapeDetection, bool);

	vtkSetMacro(Tolerance, double);
	vtkGetMacro(Tolerance, double);

	vtkSetMacro(ResultType, int);
	vtkGetMacro(ResultType, int);


	vtkSetMacro(Smoothing, bool);
	vtkGetMacro(Smoothing, bool);

	vtkSetMacro(ManualSmoothing, bool);
	vtkGetMacro(ManualSmoothing, bool);

	vtkSetMacro(SigmaDivisor, double);
	vtkGetMacro(SigmaDivisor, double);

	vtkSetMacro(TestingRadius, double);
	vtkGetMacro(TestingRadius, double);

	vtkSetMacro(Normalize, bool);
	vtkGetMacro(Normalize, bool);

	vtkSetVector3Macro(SDFRadiusFactors, double);
	vtkGetVector3Macro(SDFRadiusFactors, double);

	vtkSetVector3Macro(SDFStandardDeviations, double);
	vtkGetVector3Macro(SDFStandardDeviations, double);

	vtkSetMacro(Connectivity, unsigned int);
	vtkGetMacro(Connectivity, unsigned int);


protected:
	vtkExtractMidsurface();
	~vtkExtractMidsurface();

	int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;
	int RequestInformation(vtkInformation* vtkNotUsed(request), vtkInformationVector** inputVector, vtkInformationVector* outputVector) override;
	int FillOutputPortInformation(int, vtkInformation*) override;

private:
	vtkExtractMidsurface(const vtkExtractMidsurface&) = delete;
	void operator=(const vtkExtractMidsurface&) = delete;

	// std::vector<vtkExtractCenterLine*> centerlines;

	// GUI parameters
	char* InputArray;
	unsigned int SmoothInput;
	double Threshold;
	double IntegrationStep;
	bool AutomaticStepSize;
	double RadiusFactors[3];
	double StandardDeviations[3];
	double RadiusFactorsSDF[3];
	double StandardDeviationsSDF[3];
	int ResultType;
	bool GoldenSectionSearch;
	bool ShapeDetection;
	double Tolerance;

	// SDF parameters
	unsigned int DistanceType; // 1 - 2D, 2 - 3D
	unsigned int FieldType;	   // 1 - INSIDE, 2 - OUTSIDE, 3 - BOTH
	bool Smoothing;
	bool ManualSmoothing;
	double SDFRadiusFactors[3];
	double SDFStandardDeviations[3];
	double SigmaDivisor;
	double TestingRadius;
	bool Normalize;
	unsigned int Connectivity;
};

#endif // __vtkExtractMidsurface_h
