#ifndef __vtkComputeEigenVectorField_h
#define __vtkComputeEigenVectorField_h

#include "vtkMidsurfaceExtractorModule.h" // for export
#include <vtkPVLogger.h>
#include <vtkImageData.h>
#include <vtkGradientFilter.h>
#include <vtkDataArray.h>

#include "vtkImageAlgorithm.h"

class VTKMIDSURFACEEXTRACTOR_EXPORT vtkComputeEigenVectorField : public vtkImageAlgorithm
{
public:
	enum ESmoothInput
	{
		SMOOTH_INPUT_OFF = 1,
		SMOOTH_INPUT_GAUSSIAN = 2,
		SMOOTH_INPUT_SIGNED_DISTANCE_FIELD_SIMPLE = 3,
		SMOOTH_INPUT_SIGNED_DISTANCE_FIELD_MINKOWSI = 4
	};

	static vtkComputeEigenVectorField *New();
	vtkTypeMacro(vtkComputeEigenVectorField, vtkImageAlgorithm);

	vtkGetStringMacro(InputArray);
	vtkSetStringMacro(InputArray);

	vtkSetMacro(SmoothInput, unsigned int);
	vtkGetMacro(SmoothInput, unsigned int);

	vtkSetMacro(DistanceType, unsigned int);
	vtkGetMacro(DistanceType, unsigned int);

	vtkSetMacro(FieldType, unsigned int);
	vtkGetMacro(FieldType, unsigned int);

	vtkSetMacro(ExponentDist, double);
	vtkGetMacro(ExponentDist, double);

	vtkSetMacro(Exponent, double);
	vtkGetMacro(Exponent, double);

	vtkSetMacro(Threshold, double);
	vtkGetMacro(Threshold, double);

	vtkSetVector3Macro(RadiusFactors, double);
	vtkGetVector3Macro(RadiusFactors, double);

	vtkSetVector3Macro(StandardDeviations, double);
	vtkGetVector3Macro(StandardDeviations, double);

	vtkSetVector3Macro(RadiusFactorsSDF, double);
	vtkGetVector3Macro(RadiusFactorsSDF, double);

	vtkSetVector3Macro(StandardDeviationsSDF, double);
	vtkGetVector3Macro(StandardDeviationsSDF, double);
	
protected:
	vtkComputeEigenVectorField();
	~vtkComputeEigenVectorField();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:
    vtkComputeEigenVectorField(const vtkComputeEigenVectorField &) = delete;
    void operator=(const vtkComputeEigenVectorField &) = delete;

    void ComputeGaussianSmoothing(vtkImageData *output, vtkGradientFilter *grad);
    void ComputeSimpleSignedDistanceField(vtkImageData *output, vtkGradientFilter *grad);
    void ComputeSignedDistanceField(vtkImageData *output, vtkGradientFilter *grad);
    void ComputeVectorField(vtkDataArray *gradArr, vtkDataArray *tensArr);
    void AddRegionIDs(vtkImageData *output);
	void dfs(int x, int y, int c, std::vector<std::vector<int>> &g, std::vector<std::vector<int>> &w);

    template <class TReal>
    TReal **create_matrix(long nrow, long ncol);
    template <class TReal>
	void free_matrix(TReal **m);

	char *InputArray;
	unsigned int SmoothInput;
	unsigned int DistanceType; // 1 - MINKOWSKI, 2 - CHEBYSHEV, 3 - QUADRATIC, 4 - OCTILE, 5 - CHAMFER
	unsigned int FieldType;	   // 1 - INSIDE, 2 - OUTSIDE, 3 - BOTH
	double Exponent;	   // 1 = Manhattan distance, 2 = Euclidean distance
	double RadiusFactors[3];
	double StandardDeviations[3];
	double RadiusFactorsSDF[3];
	double StandardDeviationsSDF[3];
	double ExponentDist;
	double Threshold;
};

#endif // __vtkComputeEigenVectorField_h
