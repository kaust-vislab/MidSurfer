#ifndef __vtkExtractCenterLine_h
#define __vtkExtractCenterLine_h

#include "vtkMidSurferModule.h" // for export
#include <vtkPVLogger.h>

#include "vtkImageAlgorithm.h"
#include <vtkPoints.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkAppendPolyData.h>
#include <vtkImageInterpolator.h>
#include <vtkSmartPointer.h>

typedef std::array<double, 3> Point3;
typedef std::array<double, 3> Vector3;

namespace Midsurfacer
{
	namespace Tools
	{
		enum EResultType
		{
			RESULT_TYPE_TRIANGULATION = 1,
			RESULT_TYPE_LINE_SET = 2,
			RESULT_TYPE_POINT_SET = 3
		};

		template <class TReal>
		TReal **create_matrix(const long nrow, const long ncol);
		template <class TReal>
		void free_matrix(TReal **m);

		class CenterlineFromRegionExtractor
		{
		public:
			CenterlineFromRegionExtractor();
			~CenterlineFromRegionExtractor();

			void SetInputArray(std::string val) { this->InputArray = val; }
			void SetGoldenSectionSearch(bool val) { this->GoldenSectionSearch = val; }
			void SetResultType(int val) { this->ResultType = val; }
			void SetIntegrationStep(double val) { this->IntegrationStep = val; }
			void SetTolerance(double val) { this->Tolerance = val; }

			void ExtractCenterlineFromRegion(const Point3 &arr, vtkImageData *input, vtkPolyData *centerline);

		private:
			CenterlineFromRegionExtractor(const CenterlineFromRegionExtractor &copy_from) = delete;
			CenterlineFromRegionExtractor &operator=(const CenterlineFromRegionExtractor &copy_from) = delete;

			void ComputeNextPoint(Point3 &p, Vector3 &v, Vector3 &vold);
			void InsertNextCell(vtkCellArray *lines, const vtkIdType id1, const vtkIdType id2);
			int AppendPoints(Point3 &p, int k, vtkImageData *input, vtkPoints *points, vtkCellArray *lines);
			void ComputeEigenvector(Vector3 &v, const vtkIdType id);
			void DoGoldenSectionSearch(Point3 &p, Vector3 &v);
			template <class F>
			void FindBounds(Point3 &p, Vector3 &v, double *ax, double *cx, F f);
			template <class F>
			double golden(Point3 &p, Vector3 &v, double ax, double bx, double cx, F f);
			double GetSmoothedValue(Point3 &p);
			double GetSmoothedValueAlongLine(Point3 &p, Vector3 &v, double x);

			std::string InputArray;
			vtkDataArray *heightArr;
			vtkDataArray *smoothHeightArr;
			vtkDataArray *tensArr;
			std::vector<Point3> StartPoints;
			double **a;
			double **w;
			double CurrentDirection;
			vtkSmartPointer<vtkImageInterpolator> interp;
			vtkSmartPointer<vtkImageInterpolator> interpMask;
			bool GoldenSectionSearch;
			int ResultType;
			double IntegrationStep;
			double Tolerance;
		};
	}
}

class VTKMIDSURFER_EXPORT vtkExtractCenterLine : public vtkImageAlgorithm
{
public:
	/// @brief Matches vtkExtractMidsurface::SmoothInput, but Gaussian 3D is not applicable
	enum ESmoothingMethod
	{
		SMOOTH_INPUT_OFF = 1,
		SMOOTH_INPUT_GAUSSIAN = 2,
		SMOOTH_INPUT_SIGNED_DISTANCE_FIELD_SIMPLE = 4,
		SMOOTH_INPUT_SIGNED_DISTANCE_FIELD = 5,
	};

	static vtkExtractCenterLine *New();
	vtkTypeMacro(vtkExtractCenterLine, vtkImageAlgorithm);

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

	vtkSetMacro(Connectivity, unsigned int);
	vtkGetMacro(Connectivity, unsigned int);

protected:
	vtkExtractCenterLine();
	~vtkExtractCenterLine();

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
	int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;
	int FillOutputPortInformation(int, vtkInformation *) override;

private:
	vtkExtractCenterLine(const vtkExtractCenterLine &) = delete;
	void operator=(const vtkExtractCenterLine &) = delete;

	void ExtractCenterlineFromSlice(vtkImageData *input, vtkAppendPolyData *append);
	bool IsDiskShaped(vtkImageData *slice, int regionID);
	int SelectStartingPoints(vtkImageData *slice);

	std::vector<Point3> StartPoints;

	// GUI parameters
	char *InputArray;
	unsigned int SmoothInput;
	unsigned int DistanceType; // 1 - 2D, 2 - 3D
	unsigned int FieldType;	   // 1 - INSIDE, 2 - OUTSIDE, 3 - BOTH
	double ExponentDist;	   // 1 = Manhattan distance, 2 = Euclidean distance
	double Exponent;
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
	unsigned int Connectivity;
};

#endif // __vtkExtractCenterLine_h
