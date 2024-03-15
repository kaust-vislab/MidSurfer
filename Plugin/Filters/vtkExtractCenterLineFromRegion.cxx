#include "vtkExtractCenterLineFromRegion.h"
#include <vtkImageData.h>
#include <vtkGradientFilter.h>
#include <vtkArrayCalculator.h>
#include <vtkPassSelectedArrays.h>
#include <vtkDataArraySelection.h>
#include <vtkPassArrays.h>
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
#include <vtkExtractVOI.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include <vtkContourFilter.h>
#include <vtkPointLocator.h>
#include <vtkNew.h>
#include <vtkFeatureEdges.h>
#include <vtkCenterOfMass.h>
#include <vtkGeometryFilter.h>

#include "vtkSignedDistanceField.h"

#include <vtkXMLImageDataWriter.h>

#include <chrono>
#include <numbers>

#define SHFT2(a, b, c) \
	(a) = (b);         \
	(b) = (c);
#define SHFT3(a, b, c, d) \
	(a) = (b);            \
	(b) = (c);            \
	(c) = (d);

vtkStandardNewMacro(vtkExtractCenterLineFromRegion);

//----------------------------------------------------------------------------
vtkExtractCenterLineFromRegion::vtkExtractCenterLineFromRegion() : heightArr(nullptr), tensArr(nullptr), InputArray(nullptr)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLineFromRegion::vtkExtractCenterLineFromRegion() BEGIN");

	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	a = create_matrix<double>(3, 3);
	w = create_matrix<double>(3, 3);

	interp = vtkSmartPointer<vtkImageInterpolator>::New();
	interp->SetInterpolationModeToCubic();

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLineFromRegion::vtkExtractCenterLineFromRegion() END");
}

vtkExtractCenterLineFromRegion::~vtkExtractCenterLineFromRegion()
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLineFromRegion::~vtkExtractCenterLineFromRegion() BEGIN");

	free_matrix(w);
	free_matrix(a);

	tensArr = nullptr;
	// gradArr = nullptr;
	heightArr = nullptr;

	this->SetInputArray(nullptr);

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLineFromRegion::~vtkExtractCenterLineFromRegion() END");
}

int vtkExtractCenterLineFromRegion::FillOutputPortInformation(int, vtkInformation *info)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLineFromRegion::FillOutputPortInformation() BEGIN");

	info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLineFromRegion::FillOutputPortInformation() END");

	return 1;
}

int vtkExtractCenterLineFromRegion::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
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

int vtkExtractCenterLineFromRegion::RequestData(vtkInformation *vtkNotUsed(request),
												vtkInformationVector **inputVector,
												vtkInformationVector *outputVector)
{
	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLineFromRegion::RequestData() BEGIN");

	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and ouptut
	auto input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	auto output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkNew<vtkXMLImageDataWriter> writer;
	writer->SetInputData(input);
	writer->SetFileName("/Users/theusst/test.vti");
	writer->Write();
	
	int dims[3];
	double spacing[3];
	double bounds[6];

	input->GetDimensions(dims);
	input->GetSpacing(spacing); // not used yet, but we might want to compute the integration step size based on it
	input->GetBounds(bounds);	// not used...

	heightArr = input->GetPointData()->GetArray(this->InputArray);
	smoothHeightArr = input->GetPointData()->GetArray("smooth");
	// gradArr = input->GetPointData()->GetArray("Gradient");
	tensArr = input->GetPointData()->GetArray("T");
		
	if (this->AutomaticStepSize)
		this->IntegrationStep = std::numbers::sqrt2 * spacing[2];


	vtkNew<vtkPassSelectedArrays> pass;
	pass->SetInputData(input);
	pass->GetPointDataArraySelection()->EnableArray("smooth");
	pass->Update();

	interp->Initialize(pass->GetOutput());
	interp->Update();

	if (dims[2] == 1)
	{
		Point3 p{this->StartPoint[0], this->StartPoint[1], this->StartPoint[2]};
		vtkNew<vtkPolyData> centerline;
		this->ExtractCenterlineFromRegion(p, input, centerline);
		output->ShallowCopy(centerline);
	}
	else
	{
		vtkErrorMacro("vtkExtractCenterlineFromRegion only works on images (dim z = 1).");
	}

	vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkExtractCenterLineFromRegion::RequestData() END");

	return 1;
}

template <class TReal>
TReal **vtkExtractCenterLineFromRegion::create_matrix(const long nrow, const long ncol)
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
void vtkExtractCenterLineFromRegion::free_matrix(TReal **m)
{
	free(m[0]);
	delete[] m;
}

void vtkExtractCenterLineFromRegion::ExtractCenterlineFromRegion(const Point3 &arr, vtkImageData *slice, vtkPolyData *centerline)
{

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

	if (this->ResultType == 1) // line set
	{
		centerline->SetLines(lines);
	}
	else if (this->ResultType == 2) // point set
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

double vtkExtractCenterLineFromRegion::GetSmoothedValueAlongLine(const Point3 &p, const Vector3 &v, double x)
{
	Point3 p1;
	for (int k : {0, 1, 2})
		p1[k] = p[k] + x * v[k]; // p1 = p + x*v
	return this->GetSmoothedValue(p1);
}

double vtkExtractCenterLineFromRegion::GetSmoothedValue(const Point3 &p)
{
	double val;
	interp->Interpolate(p.data(), &val);
	return val;
}

void vtkExtractCenterLineFromRegion::FindBounds(const Point3 &p, const Vector3 &v, double *ax, double *cx, double (vtkExtractCenterLineFromRegion::*func)(const Point3 &))
{
	double dist = 0.1 * this->IntegrationStep;
	double fbx = (this->*func)(p);

	*ax = 0.0;
	*cx = 0.0;

	Point3 p1 = p;

	do
	{
		for (int k : {0, 1, 2})
			p1[k] += dist * v[k]; // p1 = p1 + dist*v
		*cx += dist;
	} while ((this->*func)(p1) > fbx);

	Point3 p2 = p;

	do
	{
		for (int k : {0, 1, 2})
			p2[k] -= dist * v[k]; // p2 = p2 - dist*v
		*ax -= dist;
	} while ((this->*func)(p2) > fbx);
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
double vtkExtractCenterLineFromRegion::golden(Point3 &p, Vector3 &v, double ax, double bx, double cx, double (vtkExtractCenterLineFromRegion::*func)(const Point3 &, const Vector3 &, double))
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

	double f1 = (this->*func)(p, v, x1);
	double f2 = (this->*func)(p, v, x2);

	while (fabs(x3 - x0) > this->Tolerance * (fabs(x1) + fabs(x2)))
	{
		if (f2 > f1)
		{
			SHFT3(x0, x1, x2, R * x1 + C * x3)
			SHFT2(f1, f2, (this->*func)(p, v, x2))
		}
		else
		{
			SHFT3(x3, x2, x1, R * x2 + C * x0)
			SHFT2(f2, f1, (this->*func)(p, v, x1))
		}
	}

	return (f1 > f2) ? x1 : x2;
}

void vtkExtractCenterLineFromRegion::DoGoldenSectionSearch(Point3 &p, Vector3 &v)
{
	Vector3 vp = {-v[1], v[0], 0.0};

	double ax, cx;
	double bx = 0.0;

	this->FindBounds(p, vp, &ax, &cx, &vtkExtractCenterLineFromRegion::GetSmoothedValue);
	double xmin = this->golden(p, vp, ax, bx, cx, &vtkExtractCenterLineFromRegion::GetSmoothedValueAlongLine);

	p[0] = p[0] + xmin * vp[0];
	p[1] = p[1] + xmin * vp[1];
	p[2] = p[2] + xmin * vp[2];
}

void vtkExtractCenterLineFromRegion::ComputeNextPoint(Point3 &p, Vector3 &v, Vector3 &vold)
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
void vtkExtractCenterLineFromRegion::InsertNextCell(vtkCellArray *lines, const vtkIdType id1, const vtkIdType id2)
{
	vtkNew<vtkLine> line;
	line->GetPointIds()->SetId(0, id1);
	line->GetPointIds()->SetId(1, id2);
	lines->InsertNextCell(line);
}

int vtkExtractCenterLineFromRegion::AppendPoints(Point3 &p, int k, vtkImageData *slice, vtkPoints *points, vtkCellArray *lines)
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

	return loop ? 0 : k + 1;
}

/// @brief Computes the eigenvector of a symmetric matrix that corresponds to the smaller eigenvalue in the plane.
/// @param v The double array the eigenvector will be stored in.
/// @param tensArr A vtkDataArray of symmetric 3x3 matrices.
/// @param k The index of the tensor in the data array.
void vtkExtractCenterLineFromRegion::ComputeEigenvector(Vector3 &v, const vtkIdType id)
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