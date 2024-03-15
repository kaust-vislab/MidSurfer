#include "vtkExtractCenterLine.h"
#include <vtkImageData.h>
#include <vtkGradientFilter.h>
#include <vtkArrayCalculator.h>
#include <vtkPassSelectedArrays.h>
#include <vtkDataArraySelection.h>
#include <vtkAppendFilter.h>
#include <vtkExtractVectorComponents.h>
#include <vtkDataSetAttributes.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkVector.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkStreamTracer.h>
#include <vtkLineSource.h>
#include <vtkLine.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include <vtkContourFilter.h>
#include <vtkPointLocator.h>
#include <vtkNew.h>
#include <vtkFeatureEdges.h>
#include <vtkCenterOfMass.h>
#include <vtkGeometryFilter.h>
#include <vtkImageMathematics.h>
#include <vtkConnectivityFilter.h>

#include "vtkSignedDistanceField.h"
#include "vtkConnectedCommponentsBinaryImage.h"

#include <chrono>
#include <numbers>

vtkStandardNewMacro(vtkExtractCenterLine);

namespace MidsurfaceExtractor
{
	namespace Tools
	{

#define SHFT2(a, b, c) \
	(a) = (b);         \
	(b) = (c);
#define SHFT3(a, b, c, d) \
	(a) = (b);            \
	(b) = (c);            \
	(c) = (d);

		CenterlineFromRegionExtractor::CenterlineFromRegionExtractor() : heightArr(nullptr), smoothHeightArr(nullptr), tensArr(nullptr), InputArray("")
		{
			a = create_matrix<double>(3, 3);
			w = create_matrix<double>(3, 3);

			interp = vtkSmartPointer<vtkImageInterpolator>::New();
			interp->SetInterpolationModeToCubic();

			interpMask = vtkSmartPointer<vtkImageInterpolator>::New();
			interpMask->SetInterpolationModeToNearest();
		}

		CenterlineFromRegionExtractor::~CenterlineFromRegionExtractor()
		{
			free_matrix(w);
			free_matrix(a);

			tensArr = nullptr;
			smoothHeightArr = nullptr;
			heightArr = nullptr;

			// this->SetInputArray(nullptr);
		}

		template <class TReal>
		TReal **CenterlineFromRegionExtractor::create_matrix(const long nrow, const long ncol)
		{
			typedef TReal *TRealPointer;
			TReal **m = new TRealPointer[nrow];

			TReal *block = (TReal *)calloc(nrow * ncol, sizeof(TReal));
			m[0] = block;
			for (int row = 1; row < nrow; ++row)
			{
				m[row] = &block[row * ncol];
			}
			return m;
		}

		/* free a TReal matrix allocated with matrix() */
		template <class TReal>
		void CenterlineFromRegionExtractor::free_matrix(TReal **m)
		{
			free(m[0]);
			delete[] m;
		}

		void CenterlineFromRegionExtractor::ExtractCenterlineFromRegion(const Point3 &arr, vtkImageData *slice, vtkPolyData *centerline)
		{

			this->heightArr = slice->GetPointData()->GetArray(this->InputArray.c_str());
			this->smoothHeightArr = slice->GetPointData()->GetArray("smooth");
			this->tensArr = slice->GetPointData()->GetArray("T");

			vtkNew<vtkPassSelectedArrays> passMask;
			passMask->SetInputData(slice);
			passMask->GetPointDataArraySelection()->EnableArray(this->InputArray.c_str());
			passMask->Update();

			passMask->GetOutput()->GetAttributes(vtkDataObject::POINT)->SetActiveScalars(this->InputArray.c_str());

			interpMask->Initialize(passMask->GetOutput());
			interpMask->Update();

			vtkNew<vtkPassSelectedArrays> pass;
			pass->SetInputData(slice);
			pass->GetPointDataArraySelection()->EnableArray("smooth");
			pass->Update();

			interp->Initialize(pass->GetOutput());
			interp->Update();

			Point3 p = arr;

			vtkLog(INFO, "Extracting line from " << p[0] << ", " << p[1] << ", " << p[2]);

			vtkNew<vtkPoints> points;
			vtkNew<vtkCellArray> lines;

			if (this->GoldenSectionSearch)
			{
				vtkIdType id = slice->FindPoint(p.data());
				Vector3 v;
				this->ComputeEigenvector(v, id);
				this->DoGoldenSectionSearch(p, v);
			}

			vtkLog(INFO, "Corrected to " << p[0] << ", " << p[1] << ", " << p[2]);

			points->InsertNextPoint(p.data());

			this->CurrentDirection = 1.0; // forward
			int k = AppendPoints(p, 1, slice, points, lines);

			if (k > 0) // if k == 0 we have a loop and are done
			{
				vtkLog(INFO, "No loop found, integrating backwards.");
				p = arr;
				this->CurrentDirection = -1.0; // backward
				AppendPoints(p, k, slice, points, lines);
			}
			else
			{
				vtkLog(INFO, "Found loop. Not integrating backwards.");
			}

			centerline->SetPoints(points);

			if ((this->ResultType == RESULT_TYPE_LINE_SET) && (centerline->GetNumberOfPoints() > 1)) // line set
			{
				centerline->SetLines(lines);
			}
			else if (this->ResultType == RESULT_TYPE_POINT_SET) // point set
			{
				vtkNew<vtkCellArray> cells;
				for (int i = 0; i < centerline->GetNumberOfPoints(); i++)
				{
					vtkIdType pid[1];
					pid[0] = i;
					cells->InsertNextCell(1, pid);
				}
				centerline->SetVerts(cells);
			}
		}

		double CenterlineFromRegionExtractor::GetSmoothedValueAlongLine(Point3 &p, Vector3 &v, double x)
		{
			Point3 p1;
			for (int k : {0, 1, 2})
				p1[k] = p[k] + x * v[k]; // p1 = p + x*v
			return this->GetSmoothedValue(p1);
		}

		double CenterlineFromRegionExtractor::GetSmoothedValue(Point3 &p)
		{
			double val;
			interp->Interpolate(p.data(), &val);
			return val;
		}

		template <class F>
		void CenterlineFromRegionExtractor::FindBounds(Point3 &p, Vector3 &v, double *ax, double *cx, F f)
		{
			double dist = 0.1 * this->IntegrationStep;
			double fbx = f(p);

			*ax = 0.0;
			*cx = 0.0;

			Point3 p1 = p;
			double val;
			do
			{
				for (int k : {0, 1, 2})
					p1[k] += dist * v[k]; // p1 = p1 + dist*v

				*cx += dist;
				interpMask->Interpolate(p1.data(), &val);
			} while ((val > 0.5) && (f(p1) > fbx));

			Point3 p2 = p;

			do
			{
				for (int k : {0, 1, 2})
					p2[k] -= dist * v[k]; // p2 = p2 - dist*v
				*ax -= dist;
				interpMask->Interpolate(p2.data(), &val);
			} while ((val > 0.5) && (f(p2) > fbx));
		}

		/// @brief Do golden section search.
		/// @param p Start point.
		/// @param v Start direction.
		/// @param ax Parameter value of left point.
		/// @param bx Parameter value of center point.
		/// @param cx Parameter value of right point.
		/// @return Parameter of maximum
		///
		/// adapted from:
		/// http://www.foo.be/docs-free/Numerical_Recipe_In_C/c10-1.pdf
		template <class F>
		double CenterlineFromRegionExtractor::golden(Point3 &p, Vector3 &v, double ax, double bx, double cx, F f)
		{
			const double R = 0.61803399;
			const double C = (1.0 - R);
			double x0, x1, x2, x3;

			x0 = ax;
			x3 = cx;
			if (fabs(cx - bx) > fabs(bx - ax))
			{
				x1 = bx;
				x2 = bx + C * (cx - bx);
			}
			else
			{
				x2 = bx;
				x1 = bx - C * (bx - ax);
			}

			double f1 = f(p, v, x1);
			double f2 = f(p, v, x2);

			while (fabs(x3 - x0) > this->Tolerance * (fabs(x1) + fabs(x2)))
			{
				if (f2 > f1)
				{
					SHFT3(x0, x1, x2, R * x1 + C * x3)
					SHFT2(f1, f2, f(p, v, x2))
				}
				else
				{
					SHFT3(x3, x2, x1, R * x2 + C * x0)
					SHFT2(f2, f1, f(p, v, x1))
				}
			}

			return (f1 > f2) ? x1 : x2;
		}

		void CenterlineFromRegionExtractor::DoGoldenSectionSearch(Point3 &p, Vector3 &v)
		{
			Vector3 vp = {-v[1], v[0], 0.0};

			double ax, cx;
			double bx = 0.0;

			using namespace std::placeholders;

			this->FindBounds(p, vp, &ax, &cx, std::bind(&CenterlineFromRegionExtractor::GetSmoothedValue, this, _1));
			double xmin = this->golden(p, vp, ax, bx, cx, std::bind(&CenterlineFromRegionExtractor::GetSmoothedValueAlongLine, this, _1, _2, _3));

			p[0] = p[0] + xmin * vp[0];
			p[1] = p[1] + xmin * vp[1];
			p[2] = p[2] + xmin * vp[2];
		}

		void CenterlineFromRegionExtractor::ComputeNextPoint(Point3 &p, Vector3 &v, Vector3 &vold)
		{

			if ((v[0] * vold[0] + v[1] * vold[1]) < 0.0)
			{
				this->CurrentDirection *= -1.0;
			}

			p[0] = p[0] + this->CurrentDirection * this->IntegrationStep * v[0];
			p[1] = p[1] + this->CurrentDirection * this->IntegrationStep * v[1];
			p[2] = p[2] + this->CurrentDirection * this->IntegrationStep * v[2];

			if (this->GoldenSectionSearch)
				this->DoGoldenSectionSearch(p, v);

			vtkMath::Assign(v, vold);
		}

		/// @brief Add connectivity for line segment.
		/// @param lines The vtkCellArray to add the connectivity to.
		/// @param id1 The id of the first point of the line segment.
		/// @param id2 The id of the second point of the line segment.
		void CenterlineFromRegionExtractor::InsertNextCell(vtkCellArray *lines, const vtkIdType id1, const vtkIdType id2)
		{
			vtkNew<vtkLine> line;
			line->GetPointIds()->SetId(0, id1);
			line->GetPointIds()->SetId(1, id2);
			lines->InsertNextCell(line);
		}

		int CenterlineFromRegionExtractor::AppendPoints(Point3 &p, int k, vtkImageData *slice, vtkPoints *points, vtkCellArray *lines)
		{
			double curval;
			bool loop = false;
			int k0 = k;

			vtkIdType id = slice->FindPoint(p.data());
			Vector3 v;
			this->ComputeEigenvector(v, id);

			Vector3 vold;
			vold = v;

			ComputeNextPoint(p, v, vold);
			id = slice->FindPoint(p.data());
			curval = this->heightArr->GetTuple1(id);

			if (curval > 0.5) // nothing to do when outside segmentation mask
			{
				points->InsertNextPoint(p[0], p[1], p[2]);
				InsertNextCell(lines, 0, k);

				this->ComputeEigenvector(v, id);
				ComputeNextPoint(p, v, vold);
				id = slice->FindPoint(p.data());
				curval = this->heightArr->GetTuple1(id);

				while ((curval > 0.5) && (loop == false) && (k < 10000)) // stop when outside of segmentation mask or found a loop
				{
					k++;

					points->InsertNextPoint(p[0], p[1], p[2]);
					InsertNextCell(lines, k - 1, k);

					if (((k - k0) > 2) &&
						(vtkMath::Distance2BetweenPoints(points->GetPoint(0), p) < 2.0 * this->IntegrationStep)) // TODO: check if point is close to starting point
					{
						InsertNextCell(lines, k, 0);
						loop = true;
					}

					this->ComputeEigenvector(v, id);
					ComputeNextPoint(p, v, vold);
					id = slice->FindPoint(p.data());
					curval = this->heightArr->GetTuple1(id);
				}
			}
			else
				return k;

			return loop ? 0 : k + 1;
		}

		/// @brief Computes the eigenvector of a symmetric matrix that corresponds to the smaller eigenvalue in the plane.
		/// @param v The double array the eigenvector will be stored in.
		/// @param tensArr A vtkDataArray of symmetric 3x3 matrices.
		/// @param k The index of the tensor in the data array.
		void CenterlineFromRegionExtractor::ComputeEigenvector(Vector3 &v, const vtkIdType id)
		{
			double e[3];
			double *tens = tensArr->GetTuple9(id);

			a[0][0] = tens[0];
			a[0][1] = tens[1];
			a[0][2] = tens[2];
			a[1][0] = tens[3];
			a[1][1] = tens[4];
			a[1][2] = tens[5];
			a[2][0] = tens[6];
			a[2][1] = tens[7];
			a[2][2] = tens[8];

			vtkMath::Jacobi(a, e, w);

			// if the z component of the smallest eigenvector is greater than zero, it points out of the plane,
			// that is not the eigenvector we want, so we take the next one
			// (the result from vtkMath::Jacobi is sorted by increasing eigenvalues)
			int num = (w[2][0] > 0.0) ? 1 : 0;

			v[0] = w[0][num];
			v[1] = w[1][num];
			v[2] = 0.0;

			vtkMath::Normalize(v.data());
		}
	}
}

//----------------------------------------------------------------------------
vtkExtractCenterLine::vtkExtractCenterLine() : InputArray(nullptr)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLine::vtkExtractCenterLine() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLine::vtkExtractCenterLine() END");
}

vtkExtractCenterLine::~vtkExtractCenterLine()
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLine::~vtkExtractCenterLine() BEGIN");

	this->SetInputArray(nullptr);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLine::~vtkExtractCenterLine() END");
}

int vtkExtractCenterLine::FillOutputPortInformation(int, vtkInformation *info)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLine::FillOutputPortInformation() BEGIN");

	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLine::FillOutputPortInformation() END");

	return 1;
}

int vtkExtractCenterLine::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *outInfo = outputVector->GetInformationObject(0);
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

	int extent[6];
	double origin[3];
	double spacing[3];

	inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent);
	inInfo->Get(vtkDataObject::SPACING(), spacing);
	inInfo->Get(vtkDataObject::ORIGIN(), origin);

	outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
	outInfo->Set(vtkDataObject::SPACING(), spacing, 3);
	outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);

	return 1;
}

int vtkExtractCenterLine::RequestData(vtkInformation *vtkNotUsed(request),
									  vtkInformationVector **inputVector,
									  vtkInformationVector *outputVector)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLine::RequestData() BEGIN");

	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	auto input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	auto output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkNew<vtkImageData> image;
	image->DeepCopy(input); // for initilaization

	int dims[3];
	double spacing[3];
	double bounds[6];

	image->GetDimensions(dims);
	image->GetSpacing(spacing); // not used yet, but we might want to compute the integration step size based on it
	image->GetBounds(bounds);	// not used...

	if (this->AutomaticStepSize)
		this->IntegrationStep = std::numbers::sqrt2 * spacing[2];

	vtkNew<vtkAppendPolyData> append;

	if (dims[2] == 1)
		ExtractCenterlineFromSlice(image, append);
	else
	{
		vtkErrorMacro("vtkExtractCenterline only works on images (dim z = 1).");
	}

	if (append->GetTotalNumberOfInputConnections() > 0)
	{
		append->Update();
		output->ShallowCopy(append->GetOutput());
	}

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLine::RequestData() END");

	return 1;
}

void vtkExtractCenterLine::ExtractCenterlineFromSlice(vtkImageData *slice, vtkAppendPolyData *append)
{
	this->StartPoints.clear();

	if (this->SmoothInput == SMOOTH_INPUT_OFF)
	{
		if (!slice->GetPointData()->GetArray("smooth")) {
			vtkLog(ERROR, "No smooth array found!");
		}

		slice->GetAttributes(vtkDataObject::POINT)->SetActiveScalars("smooth");
	}
	else if (this->SmoothInput == SMOOTH_INPUT_GAUSSIAN)
	{
		// copy the input array to double array
		vtkNew<vtkArrayCalculator> calc;
		calc->SetInputData(slice);
		calc->AddScalarArrayName(this->InputArray);
		calc->SetFunction(this->InputArray);
		calc->SetResultArrayName("smooth");
		calc->SetResultArrayType(VTK_DOUBLE);
		calc->Update();

		vtkNew<vtkImageGaussianSmooth> smooth;
		smooth->SetInputData(calc->GetOutput());
		smooth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "smooth");
		smooth->SetRadiusFactors(this->RadiusFactors);
		smooth->SetStandardDeviations(this->StandardDeviations);
		smooth->Update();

		slice->GetPointData()->AddArray(smooth->GetOutput()->GetPointData()->GetArray("smooth"));

		slice->GetAttributes(vtkDataObject::POINT)
			->SetActiveScalars("smooth");
	}
	else if (this->SmoothInput == SMOOTH_INPUT_SIGNED_DISTANCE_FIELD_SIMPLE)
	{
		// copy the input array to double array
		vtkNew<vtkArrayCalculator> calc;
		calc->SetInputData(slice);
		calc->AddScalarArrayName(this->InputArray);
		calc->SetFunction(this->InputArray);
		calc->SetResultArrayName("smooth");
		calc->SetResultArrayType(VTK_DOUBLE);
		calc->Update();

		vtkNew<vtkImageGaussianSmooth> smooth;
		smooth->SetInputData(calc->GetOutput());
		smooth->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "smooth");
		smooth->SetRadiusFactors(this->RadiusFactorsSDF);
		smooth->SetStandardDeviations(this->StandardDeviationsSDF);
		smooth->Update();

		// Create the isosurface
		vtkNew<vtkContourFilter> contour;
		contour->SetInputConnection(smooth->GetOutputPort());
		contour->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, this->InputArray);
		contour->GenerateValues(1, this->Threshold, this->Threshold);
		contour->ComputeNormalsOff();
		contour->ComputeGradientsOff();
		contour->Update();

		if (contour->GetOutput()->GetNumberOfPoints() > 0)
		{
			// Initialize the locator
			vtkNew<vtkPointLocator> pointTree;
			pointTree->SetDataSet(contour->GetOutput());
			pointTree->BuildLocator();

			auto inArr = smooth->GetOutput()->GetPointData()->GetArray("smooth");
			double p[3];  // for voxel location
			double pc[3]; // for clostest point of contour

			double min = VTK_DOUBLE_MAX;
			double max = VTK_DOUBLE_MIN;

			for (auto i = 0; i < inArr->GetNumberOfTuples(); i++)
			{
				slice->GetPoint(i, p);
				auto val = inArr->GetTuple1(i);

				contour->GetOutput()->GetPoint(pointTree->FindClosestPoint(p), pc);
				double udist = std::sqrt(vtkMath::Distance2BetweenPoints(p, pc));
				udist = val < this->Threshold ? -1.0 * udist : udist;
				inArr->SetTuple1(i, udist);

				if (udist < min)
					min = udist;

				if (udist > max)
					max = udist;
			}

			for (auto i = 0; i < inArr->GetNumberOfTuples(); i++)
			{
				auto val = inArr->GetTuple1(i);

				if (val >= 0.0)
					inArr->SetTuple1(i, max * std::pow(val / max, this->ExponentDist));
				else
					inArr->SetTuple1(i, min * std::pow(val / min, this->ExponentDist));
			}
		}

		slice->GetPointData()->AddArray(smooth->GetOutput()->GetPointData()->GetArray("smooth"));

		slice->GetAttributes(vtkDataObject::POINT)
			->SetActiveScalars("smooth");
	}
	else if (this->SmoothInput == SMOOTH_INPUT_SIGNED_DISTANCE_FIELD)
	{
		vtkNew<vtkSignedDistanceField> dist;
		dist->SetInputData(slice);
		dist->SetInputArray(this->InputArray);
		dist->SetDistanceType(this->DistanceType);
		dist->SetFieldType(this->FieldType);
		dist->Update();

		slice->GetPointData()->AddArray(dist->GetOutput()->GetPointData()->GetArray("smooth"));

		slice->GetAttributes(vtkDataObject::POINT)
			->SetActiveScalars("smooth");
	}
	else
	{
		vtkErrorMacro("Unknown smoothing method");
		return;
	}

	if (this->SelectStartingPoints(slice) == -1 || this->StartPoints.size() == 0)
	{
		vtkLog(INFO, "No starting points found. Aborting.");
		return;
	}

	// Compute gradients
	vtkNew<vtkGradientFilter> grad;
	grad->SetInputData(slice);
	grad->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "smooth");
	grad->SetResultArrayName("Gradient");
	grad->Update();

	// Compute tensors
	vtkNew<vtkGradientFilter> tens;
	tens->SetInputConnection(grad->GetOutputPort());
	tens->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS, "Gradient");
	tens->SetResultArrayName("T");
	tens->Update();

	MidsurfaceExtractor::Tools::CenterlineFromRegionExtractor extract;

	for (Point3 p : this->StartPoints)
	{
		extract.SetInputArray(std::string(this->InputArray));
		extract.SetGoldenSectionSearch(this->GoldenSectionSearch);
		extract.SetIntegrationStep(this->IntegrationStep);
		extract.SetResultType(this->ResultType);
		extract.SetTolerance(this->Tolerance);

		vtkNew<vtkPolyData> centerline;
		extract.ExtractCenterlineFromRegion(p, vtkImageData::SafeDownCast(tens->GetOutput()), centerline);
		append->AddInputData(centerline);
	}

	append->Update();
}

int vtkExtractCenterLine::SelectStartingPoints(vtkImageData *slice)
{
	vtkNew<vtkConnectedCommponentsBinaryImage> conn;
	conn->SetInputData(slice);
	conn->SetInputArray(this->InputArray);
	conn->SetConnectivity(this->Connectivity);
	conn->Update();

	auto numberOfRegions = conn->GetNumberOfConnectedRegions();
	vtkLog(INFO, "Found " << numberOfRegions << " regions.");

	if (numberOfRegions == 0)
	{
		vtkLog(INFO, "No regions found. Aborting.");
		return -1;
	}

	auto remove = std::vector<int>();

	std::unordered_map<int, double> maxValues;
	std::unordered_map<int, Point3> startingPoints;

	if (this->ShapeDetection)
	{
		auto mask = conn->GetOutput();
		for (int k = 1; k <= numberOfRegions; k++)
		{
			if (this->IsDiskShaped(mask, k))
			{
				vtkLog(INFO, "Detected disk shape. Discarding.")
					remove.push_back(k);
			}
		}
	}

	// go through the resulting image and find the maximum value in each region
	for (int i = 0; i < slice->GetNumberOfPoints(); i++)
	{
		auto regionId = conn->GetOutput()->GetPointData()->GetArray("RegionID")->GetTuple1(i);
		// check if the region is in the list of regions to remove
		if (std::find(remove.begin(), remove.end(), regionId) != remove.end() || regionId == 0)
			continue;

		auto value = slice->GetPointData()->GetArray("smooth")->GetTuple1(i);
		auto point = slice->GetPoint(i);

		if (value > maxValues[regionId])
		{
			startingPoints[regionId] = {point[0], point[1], point[2]};
			maxValues[regionId] = value;
		}
	}

	for (auto &[k, v] : startingPoints)
	{
		this->StartPoints.push_back(v);
	}

	vtkLog(INFO, "Selected " << this->StartPoints.size() << " starting points.");
	for (auto point : this->StartPoints)
	{
		vtkLog(INFO, "x: " << point[0] << ", y: " << point[1] << ", z: " << point[2] << ".");
	}

	return 0;
}

bool vtkExtractCenterLine::IsDiskShaped(vtkImageData *mask, int regionID)
{
	mask->GetPointData()->SetActiveScalars("RegionID");

	vtkLog(INFO, "Checking region " << regionID << ".");
	double lowThresh = regionID - 0.5;
	double highThresh = regionID + 0.5;

	vtkNew<vtkThreshold> threshReg;
	threshReg->SetInputData(mask);
	threshReg->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionID");
	threshReg->SetLowerThreshold(static_cast<double>(lowThresh));
	threshReg->SetUpperThreshold(static_cast<double>(highThresh));
	threshReg->Update();

	// if the region is empty, we can skip it
	if (threshReg->GetOutput()->GetNumberOfPoints() == 0)
	{
		vtkLog(INFO, "Region is empty. Skipping.");
		return false;
	}

	vtkNew<vtkContourFilter> contourFilter;
	contourFilter->SetInputConnection(threshReg->GetOutputPort());
	contourFilter->SetValue(0, 0.5);
	contourFilter->Update();

	vtkNew<vtkConnectivityFilter> connectivityFilter;
	connectivityFilter->SetInputConnection(contourFilter->GetOutputPort());
	connectivityFilter->SetExtractionModeToAllRegions();
	connectivityFilter->Update();
	int numberOfContours = connectivityFilter->GetNumberOfExtractedRegions();
	vtkLog(INFO, "Found " << numberOfContours << " contours.");

	vtkNew<vtkGeometryFilter> surface;
	surface->SetInputData(threshReg->GetOutput());
	surface->Update();

	vtkNew<vtkFeatureEdges> boundary;
	boundary->SetInputConnection(surface->GetOutputPort());
	boundary->BoundaryEdgesOn();
	boundary->FeatureEdgesOff();
	boundary->ManifoldEdgesOff();
	boundary->NonManifoldEdgesOff();
	boundary->ColoringOff();
	boundary->Update();

	vtkNew<vtkCenterOfMass> centerOfMass;
	centerOfMass->SetInputConnection(boundary->GetOutputPort());
	centerOfMass->SetUseScalarsAsWeights(false);
	centerOfMass->Update();

	double *center = centerOfMass->GetCenter();

	double min = VTK_DOUBLE_MAX;
	double max = VTK_DOUBLE_MIN;

	for (vtkIdType id = 0; id < boundary->GetOutput()->GetNumberOfPoints(); id++)
	{
		double *p = boundary->GetOutput()->GetPoint(id);
		const double dist = std::sqrt(vtkMath::Distance2BetweenPoints(center, p));

		if (dist < min)
			min = dist;

		if (dist > max)
			max = dist;
	}

	return ((max - min) < 0.25 * max);
}