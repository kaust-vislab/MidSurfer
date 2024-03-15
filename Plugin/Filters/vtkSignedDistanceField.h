/**
 * https://github.com/fau-vc/optimal-frames
 */

#ifndef __vtkSignedDistanceField_h
#define __vtkSignedDistanceField_h

#include "vtkMidsurfaceExtractorModule.h" // for export

#include <vtkBitArray.h>
#include <vtkDataArrayMeta.h>
#include <vtkDoubleArray.h>
#include <vtkImageAlgorithm.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkSetGet.h>
#include <vtkImageAppend.h>

class VTKMIDSURFACEEXTRACTOR_EXPORT vtkSignedDistanceField
	: public vtkImageAlgorithm {
public:
	enum EDistanceType {
		_2D = 1,
		_3D = 2
	};
	enum EFieldType { INSIDE = 1, OUTSIDE = 2, BOTH = 3 };

	static vtkSignedDistanceField* New();
	vtkTypeMacro(vtkSignedDistanceField, vtkImageAlgorithm);

	vtkGetStringMacro(InputArray);
	vtkSetStringMacro(InputArray);

	vtkGetStringMacro(SDFResultName);
	vtkSetStringMacro(SDFResultName);

	vtkGetStringMacro(FinalResultName);
	vtkSetStringMacro(FinalResultName);

	vtkSetMacro(DistanceType, unsigned int);
	vtkGetMacro(DistanceType, unsigned int);

	vtkSetMacro(FieldType, unsigned int);
	vtkGetMacro(FieldType, unsigned int);

	vtkSetMacro(Smoothing, bool);
	vtkGetMacro(Smoothing, bool);

	vtkSetMacro(ManualSmoothing, bool);
	vtkGetMacro(ManualSmoothing, bool);

	vtkSetVector3Macro(RadiusFactors, double);
	vtkGetVector3Macro(RadiusFactors, double);

	vtkSetVector3Macro(StandardDeviations, double);
	vtkGetVector3Macro(StandardDeviations, double);

	vtkSetMacro(SigmaDivisor, double);
	vtkGetMacro(SigmaDivisor, double);

	vtkSetMacro(TestingRadius, double);
	vtkGetMacro(TestingRadius, double);

	vtkSetMacro(Normalize, bool);
	vtkGetMacro(Normalize, bool);

protected:
	vtkSignedDistanceField();
	~vtkSignedDistanceField();

	int RequestData(vtkInformation*, vtkInformationVector**,
		vtkInformationVector*) override;
	int RequestInformation(vtkInformation* vtkNotUsed(request),
		vtkInformationVector** inputVector,
		vtkInformationVector* outputVector) override;

private:
	vtkSignedDistanceField(const vtkSignedDistanceField&) = delete;
	void operator=(const vtkSignedDistanceField&) = delete;

	char* InputArray;
	char* SDFResultName;
	char* FinalResultName;
	int width, height, depth;
	int* dims;
	void ComputeSignedDistanceField(vtkDoubleArray* dist, vtkImageData* originalImage, vtkImageAppend* append=nullptr);
	void ComputeStandardDeviations(vtkDoubleArray* dist);
	void ComputeRadiusFactors();
	void GaussianSmooth(vtkDoubleArray* smooth_image, vtkDoubleArray* dist);

	unsigned int DistanceType; // 1 - _2D, 2 - _3D
	unsigned int FieldType;    // 1 - INSIDE, 2 - OUTSIDE, 3 - BOTH
	bool Smoothing;
	bool ManualSmoothing;
	double RadiusFactors[3];
	double StandardDeviations[3];
	double SigmaDivisor;
	double TestingRadius;
	bool Normalize;
};

#endif // __vtkSignedDistanceField_h
