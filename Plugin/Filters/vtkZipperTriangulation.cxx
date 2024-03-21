#include "vtkZipperTriangulation.h"
#include <vtkPolyData.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkMath.h>
#include <fstream> // Include the header for file I/O
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>


vtkStandardNewMacro(vtkZipperTriangulation);

vtkZipperTriangulation::vtkZipperTriangulation()
{
    vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkZipperTriangulation::vtkZipperTriangulation() BEGIN");

    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

    vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), ":vtkZipperTriangulation:vtkZipperTriangulation() END");
}

vtkZipperTriangulation::~vtkZipperTriangulation()
{
    vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), "vtkZipperTriangulation::~vtkZipperTriangulation() BEGIN");

    vtkVLog(PARAVIEW_LOG_PLUGIN_VERBOSITY(), ":vtkZipperTriangulation:~vtkZipperTriangulation() END");
}

int vtkZipperTriangulation::RequestData(vtkInformation* vtkNotUsed(request),
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{

    // Open the text file for writing
    std::ofstream logfile("Dklogfile2.txt");
    logfile << "11111111111 ...\n"; 

    // Get the info objects
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    // Get the input and output
    vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Perform the shallow copy
    output->ShallowCopy(input);

    // Create a scalar array to store the slice numbers for points
// Include the necessary header


    //------------------------------------------------------------------------------------------------------------------------------------------------
// Create a scalar array to store the slice numbers for points
    vtkSmartPointer<vtkDoubleArray> v_sliceno = vtkSmartPointer<vtkDoubleArray>::New();
    v_sliceno->SetName("PointSliceNumber");

    // Calculate and assign the slice numbers for points
    vtkIdType numPoints = output->GetNumberOfPoints();
    vtkIdType numVertices = output->GetNumberOfVerts();

    vtkIdType numEdges = output->GetNumberOfLines();
    double prevZ = output->GetPoint(0)[2];
    int integral = 1;
    for (vtkIdType i = 0; i < numPoints; ++i) {
        double z = output->GetPoint(i)[2];
        if (z > prevZ) {
            integral++;
        }
        v_sliceno->InsertNextValue(integral);
        prevZ = z;
    }
    int slice_begin = 1;
    int slice_end= v_sliceno->GetValue(numPoints-1);
    //------------------------------------------------------------------------------------------------------------------------------------------------
    //vtkIdType numEdges = 0;
    vtkCellArray* lines = output->GetLines();
    logfile << "\n num of points  = " << numPoints;
    logfile << "\n num of Vertices = " << numVertices;
    logfile << "\n num of edges   = " << numEdges; 
//------------------------------------------------------------------------------------------------------------------------------
    // Create a scalar array to store the slice numbers for edges
    vtkSmartPointer<vtkIntArray> e_sliceno = vtkSmartPointer<vtkIntArray>::New();
    e_sliceno->SetName("EdgeSliceNumber");  
    //set e_sliceno for all edges //copy from vertices 
    e_sliceno->SetNumberOfValues(numEdges);
    for (vtkIdType e = 0; e < numEdges; ++e) { e_sliceno->SetValue(e, -1); }


    //for (vtkIdType e = 0; e < numEdges; ++e) {
    //    vtkSmartPointer<vtkIdList> edgePts = vtkSmartPointer<vtkIdList>::New();
    //    lines->GetCell(e, edgePts); 
    //        vtkIdType v1 = edgePts->GetId(0);
    //        vtkIdType v2 = edgePts->GetId(1);
    //        // You can access their coordinates using output->GetPoint(pointId) 
    //        double v_slice = v_sliceno->GetValue(v1);
    //        if (v_slice == v_sliceno->GetValue(v2)) {
    //            e_sliceno->SetValue(e, v_slice);
    //        }
    //        else {
    //            e_sliceno->SetValue(e, -1);
    //        }
    //     
    //}

     

    //------------------------------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------------------------------ 
//----------------------------------------------------------------------------------------

           /* for (vtkIdType i = 0; i < numPoints; ++i) {
                logfile << "\n v_sliceno->GetValue(i) = " << v_sliceno->GetValue(i);
            }*/
    /*vtkIdType npts;
    vtkIdType* pts;*/

    //// Traverse the lines
    //for (lines->InitTraversal(); lines->GetNextCell(npts, pts);) {
    //    // Ensure each line has at least two points (i.e., forms an edge)
    //    if (npts >= 2) {
    //        vtkIdType v1 = pts[0];
    //        vtkIdType v2 = pts[1];

    //        // Do something with v1 and v2 (e.g., store them, process them)
    //    }
    //}

            for (vtkIdType i = 0; i < numEdges; ++i) {


                //vtkSmartPointer<vtkIdList> linePts = vtkSmartPointer<vtkIdList>::New();
                //lines->GetNextCell(linePts); // Get the next line composed of vertices
                //vtkIdType v1 = linePts->GetId(0);
                //vtkIdType v2 = linePts->GetId(1);

                /*lines->GetCell(i, edgePts); 
                 vtkIdType v1 = edgePts->GetId(0);
                   vtkIdType v2 = edgePts->GetId(1);*/

                   logfile << "\n  e_sliceno[" << i << "] = " << e_sliceno->GetValue(i);
                  // logfile << "\t ===>  v1 =" << v1 << "] v2 =" << v2;
                } 

            


//----------------------------------------------------------------------------------------


            return 1;
             

        // Close the text file
        logfile.close();
     
     


    return 1;
}

int vtkZipperTriangulation::FillInputPortInformation(int port, vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}
