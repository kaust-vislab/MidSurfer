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

#include <vector>
#include <vtkTriangle.h>
#include <vtkMath.h>
#include <vtkQuad.h> 


#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>

#include <cmath>
#include <vtkPoints.h>
#include <vtkCleanPolyData.h>



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


//------------------------------------------------------------------

void repair_colocate_verticesdd(vtkSmartPointer<vtkPolyData> mesh, double colocate_epsilon) {
    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputData(mesh);
    cleaner->PointMergingOn();
    cleaner->SetTolerance(colocate_epsilon);
    cleaner->Update();
    mesh->DeepCopy(cleaner->GetOutput());
}
//------------------------------------------------
void repair_colocate_vertices(vtkPolyData* mesh, double colocate_epsilon) {
    vtkPoints* points = mesh->GetPoints();
    if (!points) {
        return;
    }

    // Initialize a map to store old vertex indices to new vertex indices
    std::map<vtkIdType, vtkIdType> old2new;

    // Initialize the number of new vertices
    vtkIdType nb_new_vertices = 0;

    // Loop through each vertex to handle colocation
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
        double point[3];
        points->GetPoint(i, point);

        // Check for colocation based on epsilon
        bool colocated = false;
        for (const auto& entry : old2new) {
            vtkIdType old_index = entry.first;
            double* old_point = points->GetPoint(old_index);
            if (vtkMath::Distance2BetweenPoints(point, old_point) <= colocate_epsilon * colocate_epsilon) {
                old2new[i] = old_index; // Map current index to existing index
                colocated = true;
                break;
            }
        }
        if (!colocated) {
            old2new[i] = nb_new_vertices; // Map current index to a new index
            nb_new_vertices++;
        }
    }

    // Check if no colocation occurred
    if (nb_new_vertices == mesh->GetNumberOfPoints()) {
        return;
    }

    // Replace vertex indices for cells
    vtkCellArray* polys = mesh->GetPolys();
    if (polys) {
        vtkIdTypeArray* cellArray = polys->GetData();
        vtkIdType* cellIndices = cellArray->GetPointer(0);
        vtkIdType numberOfCells = cellArray->GetNumberOfTuples() / (cellArray->GetNumberOfComponents() + 1);

        for (vtkIdType i = 0; i < numberOfCells; ++i) {
            vtkIdType numberOfCellPoints = cellIndices[i * (cellArray->GetNumberOfComponents() + 1)];
            for (vtkIdType j = 0; j < numberOfCellPoints; ++j) {
                vtkIdType oldIndex = cellIndices[i * (cellArray->GetNumberOfComponents() + 1) + 1 + j];
                cellIndices[i * (cellArray->GetNumberOfComponents() + 1) + 1 + j] = old2new[oldIndex];
            }
        }
    }

    // Remove the duplicated vertices
    vtkSmartPointer<vtkIdList> deleteVerts = vtkSmartPointer<vtkIdList>::New();
    for (const auto& entry : old2new) {
        if (entry.first != entry.second) {
            deleteVerts->InsertUniqueId(entry.first);
        }
    }
    for (vtkIdType i = 0; i < deleteVerts->GetNumberOfIds(); ++i) {
        mesh->DeletePoint(deleteVerts->GetId(i));
    }
    mesh->RemoveDeletedCells();
}


//--------------------------------------------------------------------------------------------------------------------------------
vtkIdType NearestVertex2v(vtkPolyData* mesh, vtkIdType v1_e1, vtkIdType e2) {
    vtkPoints* points = mesh->GetPoints();
    if (!points) {
        return 0.0;
    }

    // Get vertex IDs of the two edges 
    vtkIdType v1_e2 = mesh->GetCell(e2)->GetPointId(0);
    vtkIdType v2_e2 = mesh->GetCell(e2)->GetPointId(1);

    // Calculate midpoint of the two edges
    double p1[3], p2[3], p3[3];
    for (int i = 0; i < 3; ++i) {
        p1[i] = points->GetPoint(v1_e1)[i];
        p2[i] = points->GetPoint(v1_e2)[i];
        p3[i] = points->GetPoint(v2_e2)[i];
    }


    // Calculate distance between midpoints
    double d1 = 0;
    double d2 = 0;
    for (int i = 0; i < 3; ++i) {
        d1 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
        d2 += (p1[i] - p3[i]) * (p1[i] - p3[i]);
    }
    if (d1 < d2)
    {
        return v1_e2;
    }
    else {
        return v2_e2;

    }
}

//--------------------------------------------------------------------------------------------------------------------------------
vtkIdType NearestVertex2Mid(vtkPolyData* mesh, vtkIdType v1_e1, vtkIdType v2_e1, vtkIdType e2) {
    vtkPoints* points = mesh->GetPoints();
    if (!points) {
        return 0.0;
    }

    // Get vertex IDs of the two edges 
    vtkIdType v1_e2 = mesh->GetCell(e2)->GetPointId(0);
    vtkIdType v2_e2 = mesh->GetCell(e2)->GetPointId(1);

    // Calculate midpoint of the two edges
    double p1[3], p2[3], p3[3];
    for (int i = 0; i < 3; ++i) {
        p1[i] = (points->GetPoint(v1_e1)[i] + points->GetPoint(v2_e1)[i]) / 2.0;
        p2[i] = points->GetPoint(v1_e2)[i];
        p3[i] = points->GetPoint(v2_e2)[i];
    }


    // Calculate distance between midpoints
    double d1 = 0;
    double d2 = 0;
    for (int i = 0; i < 3; ++i) {
        d1 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
        d2 += (p1[i] - p3[i]) * (p1[i] - p3[i]);
    }
    if (d1 < d2)
    {
        return v1_e2;
    }
    else {
        return v2_e2;

    }
}

//--------------------------------------------------------------------------------------------------------------------------------
double dist_bw_2es(vtkPolyData* mesh, vtkIdType e1, vtkIdType e2) {
    vtkPoints* points = mesh->GetPoints();
    if (!points) {
        return 0.0;
    }

    // Get vertex IDs of the two edges
    vtkIdType v1_e1 = mesh->GetCell(e1)->GetPointId(0);
    vtkIdType v2_e1 = mesh->GetCell(e1)->GetPointId(1);
    vtkIdType v1_e2 = mesh->GetCell(e2)->GetPointId(0);
    vtkIdType v2_e2 = mesh->GetCell(e2)->GetPointId(1);

    // Calculate midpoint of the two edges
    double p1[3], p2[3];
    for (int i = 0; i < 3; ++i) {
        p1[i] = (points->GetPoint(v1_e1)[i] + points->GetPoint(v2_e1)[i]) / 2.0;
        p2[i] = (points->GetPoint(v1_e2)[i] + points->GetPoint(v2_e2)[i]) / 2.0;
    }

    // Calculate distance between midpoints
    double d1 = 0;
    for (int i = 0; i < 3; ++i) {
        d1 += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return sqrt(d1);
}
//...........................................................................................................................  .......................  .......................  .......................  
bool Create1Triangle(vtkPolyData* mesh, vtkSmartPointer<vtkCellArray> cells, vtkIdType e1, vtkIdType e2, vtkIdType e2btm, double d) {
    //we can use d for avoiding triangles with longer edges 
    /*if (dist_bw_2es(mesh, e1, e2) > d)
        return false;*/
    vtkPoints* points = mesh->GetPoints();
    if (!points) {
        return false;
    }
    vtkIdType v1 = mesh->GetCell(e1)->GetPointId(0);
    vtkIdType v2c = mesh->GetCell(e1)->GetPointId(1);

    //v2c is supposed to be common vertex beween e1 and e2btm 
    vtkIdType v1n = mesh->GetCell(e2btm)->GetPointId(0);
    vtkIdType v2n = mesh->GetCell(e2btm)->GetPointId(1);
    if (v1n==v1||v2n==v1)
    {
        vtkIdType v1t = v2c;
        v2c = v1;
        v1 = v1t;//swaping 
    }

    vtkIdType v3 = NearestVertex2v(mesh, v1, e2);
    vtkIdType triangle1[3] = { v1, v2c, v3 };
    cells->InsertNextCell(3, triangle1);
    return true;
}
//===================================================================
void Cross(const double v1[3], const double v2[3], double result[3]) {
result[0] = v1[1] * v2[2] - v1[2] * v2[1];
result[1] = v1[2] * v2[0] - v1[0] * v2[2];
result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void Subtract(const double v1[3], const double v2[3], double result[3]) {
    for (int i = 0; i < 3; ++i) {
        result[i] = v1[i] - v2[i];
    }
}

double Dot(const double v1[3], const double v2[3]) {
    double dotProduct = 0.0;
    for (int i = 0; i < 3; ++i) {
        dotProduct += v1[i] * v2[i];
    }
    return dotProduct;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Create2TrianglesDD(vtkPolyData* mesh, vtkSmartPointer<vtkCellArray> cells, vtkIdType e1, vtkIdType e2, double d) {
    // Retrieve points from the mesh
    vtkPoints* points = mesh->GetPoints();
    if (!points) {
        return false;
    }

    // Get the indices of the vertices connected by edges e1 and e2
    vtkIdType v1 = mesh->GetCell(e1)->GetPointId(0);
    vtkIdType v2 = mesh->GetCell(e1)->GetPointId(1);
    vtkIdType v3 = mesh->GetCell(e2)->GetPointId(0);
    vtkIdType v4 = mesh->GetCell(e2)->GetPointId(1);

    // Calculate vectors for each edge
    double p1[3], p2[3], p3[3], p4[3];
    points->GetPoint(v1, p1);
    points->GetPoint(v2, p2);
    points->GetPoint(v3, p3);
    points->GetPoint(v4, p4);

    // Calculate the difference vectors
    double s1[3], s2[3];
    Subtract(p2, p1, s1);
    Subtract(p3, p2, s2);

    // Calculate the dot products
    double dotProduct1 = Dot(s1, s1);
    double dotProduct2 = Dot(s2, s2);


    vtkIdType triangle1[3], triangle2[3];

    // Adjust the face orientations if necessary
    if (dotProduct1 > 0 && dotProduct2 > 0) { // Both cross products point outward
        // Define triangles with counterclockwise winding order
        triangle1[0] = v1; triangle1[1] = v2; triangle1[2] = v3;
        triangle2[0] = v4; triangle2[1] = v3; triangle2[2] = v2;
    }
    else if (dotProduct1 < 0 && dotProduct2 < 0) { // Both cross products point inward
     // Define triangles with counterclockwise winding order
        triangle1[0] = v2; triangle1[1] = v1; triangle1[2] = v4;
        triangle2[0] = v1; triangle2[1] = v3; triangle2[2] = v4;
    }
    else { // Adjust one triangle to ensure consistent face orientations
        if (dotProduct1 < 0) { // Cross product 1 points inward
            // Define triangles with counterclockwise winding order
            triangle1[0] = v1; triangle1[1] = v2; triangle1[2] = v3;
            triangle2[0] = v4; triangle2[1] = v3; triangle2[2] = v2;
        }
        else { // Cross product 2 points inward
         // Define triangles with counterclockwise winding order
            triangle1[0] = v2; triangle1[1] = v1; triangle1[2] = v4;
            triangle2[0] = v1; triangle2[1] = v3; triangle2[2] = v4;
        }
    }

    // Insert triangles into the cell array
    cells->InsertNextCell(3, triangle1);
    cells->InsertNextCell(3, triangle2);

    return true;
}

//...........................................................................................................................  
bool Create2Triangles(vtkPolyData* mesh, vtkSmartPointer<vtkCellArray> cells, vtkIdType e1, vtkIdType e2, double d) {
    //we can use d for avoiding triangles with longer edges 
    /*if (dist_bw_2es(mesh, e1, e2) > d)
        return false;*/
    vtkPoints* points = mesh->GetPoints();
    if (!points) {
        return false;
    }
    vtkIdType v1 = mesh->GetCell(e1)->GetPointId(0);
    vtkIdType v2 = mesh->GetCell(e1)->GetPointId(1);
    vtkIdType v3 = mesh->GetCell(e2)->GetPointId(0);
    vtkIdType v4 = mesh->GetCell(e2)->GetPointId(1);


    // Calculate vectors for each edge
    double p1[3], p2[3], p3[3], p4[3], dir1[3], dir2[3];

    for (int i = 0; i < 3; ++i) {
        p1[i] = points->GetPoint(v1)[i];
        p2[i] = points->GetPoint(v2)[i];
        p3[i] = points->GetPoint(v3)[i];
        p4[i] = points->GetPoint(v4)[i];
    }

   /* points->GetPoint(v1, p1);
    points->GetPoint(v2, p2);
    points->GetPoint(v3, p3);
    points->GetPoint(v4, p4);*/

    // Calculate direction vectors
    for (int i = 0; i < 3; ++i) {
        dir1[i] = p2[i] - p1[i];
        dir2[i] = p4[i] - p3[i];
    } 
     
    // Check the direction of edges e1 and e2
    double dotProduct = vtkMath::Dot(dir1, dir2);
    //vtkSmartPointer<vtkCellArray> cells = mesh->GetPolys();
    if (dotProduct > 0) { // Case 1 and Case 2
         
        double p1p2[3], p4p2[3];
        for (int i = 0; i < 3; ++i) {
            p1p2[i] = p1[i] - p2[i];
            p4p2[i] = p4[i] - p2[i];
        }
        double p2p1[3], p3p1[3];
        for (int i = 0; i < 3; ++i) {
            p2p1[i] = p2[i] - p1[i];
            p3p1[i] = p3[i] - p1[i];
        }


        double v1v2v4 = vtkMath::AngleBetweenVectors(p1p2, p4p2);//angle at v2 
        //This calculates the angle at vertex v2 using the vectors from v1 to v2 and from v2 to v4.
        double v2v1v3 = vtkMath::AngleBetweenVectors(p2p1, p3p1);//angle at v1
        if (v1v2v4 > v2v1v3) { // Case 1
            /*if (p1[0] <p2[0])
            {*/
                vtkIdType triangle1[3] = { v1, v2, v3 };  // Counterclockwise order
                vtkIdType triangle2[3] = { v4, v3, v2 };  // Counterclockwise order
                cells->InsertNextCell(3, triangle1);
                cells->InsertNextCell(3, triangle2);
            //}
            //else
            //{ 
            //    vtkIdType triangle1[3] = { v2, v1, v3 };  // Counterclockwise order
            //    vtkIdType triangle2[3] = { v3, v4, v2 };  // Counterclockwise order
            //    cells->InsertNextCell(3, triangle1);
            //    cells->InsertNextCell(3, triangle2);
            //}
           // vtkIdType triangle2[3] = { v4 ,v2, v3};  // dd
        }
        else { // Case 2 
            vtkIdType triangle1[3] = { v2, v1, v4 };  // Counterclockwise order
            vtkIdType triangle2[3] = { v1, v3, v4 };  // Counterclockwise order
            cells->InsertNextCell(3, triangle1);
            cells->InsertNextCell(3, triangle2);
        }
    }
    else { // Case 3 and Case 4
        // opposit direction so dignoal may either be v2v4 OR v1v3  
        double p1p2[3], p3p2[3];
        for (int i = 0; i < 3; ++i) {
            p1p2[i] = p1[i] - p2[i];
            p3p2[i] = p3[i] - p2[i];
        }
        double p2p1[3], p4p1[3];
        for (int i = 0; i < 3; ++i) {
            p2p1[i] = p2[i] - p1[i];
            p4p1[i] = p4[i] - p1[i];
        }


        double v1v2v3 = vtkMath::AngleBetweenVectors(p1p2, p3p2);//angle at v2
        double v2v1v4 = vtkMath::AngleBetweenVectors(p2p1,p4p1);//angle at v1
        if (v1v2v3 > v2v1v4) { // Case 3
            vtkIdType triangle1[3] = { v2, v1, v4 };  // Counterclockwise order
            vtkIdType triangle2[3] = { v3, v2, v4 };  // Counterclockwise order
            cells->InsertNextCell(3, triangle1);
            cells->InsertNextCell(3, triangle2);
        }
        else { // Case 4
            vtkIdType triangle1[3] = { v2, v1, v3 };  // Counterclockwise order
            vtkIdType triangle2[3] = { v1, v4, v3 };  // Counterclockwise order
            cells->InsertNextCell(3, triangle1);
            cells->InsertNextCell(3, triangle2);
        }
    }

    return true;
}

//    if (dotProduct > 0) { 
//        double p1p2[3], p4p2[3];
//        for (int i = 0; i < 3; ++i) {
//            p1p2[i] = p1[i] - p2[i];
//            p4p2[i] = p4[i] - p2[i];
//        }
//        double p2p1[3], p3p1[3];
//        for (int i = 0; i < 3; ++i) {
//            p2p1[i] = p2[i] - p1[i];
//            p3p1[i] = p3[i] - p1[i];
//        }
//
//
//        double v1v2v4 = vtkMath::AngleBetweenVectors(p1p2, p4p2);//angle at v2 
//        //This calculates the angle at vertex v2 using the vectors from v1 to v2 and from v2 to v4.
//        double v2v1v3 = vtkMath::AngleBetweenVectors(p2p1, p3p1);//angle at v1
//        if (v1v2v4 > v2v1v3) {  
//            vtkIdType triangle1[3] = { v1, v2, v3 };
//            vtkIdType triangle2[3] = { v4, v3, v2 };
//            cells->InsertNextCell(3,triangle1);
//            cells->InsertNextCell(3, triangle2);
//
//        }
//        else {
//            vtkIdType triangle1[3] = { v2, v1, v4 };
//            vtkIdType triangle2[3] = { v1, v3, v4 }; 
//            cells->InsertNextCell(3, triangle1);
//            cells->InsertNextCell(3, triangle2);
//
//            // Set points and triangles to the polyData 
//
//        }
//
//
//    }
//    else {// opposit direction so dignoal may either be v2v4 OR v1v3 
//        double p1p2[3], p3p2[3];
//        for (int i = 0; i < 3; ++i) {
//            p1p2[i] = p1[i] - p2[i];
//            p3p2[i] = p3[i] - p2[i];
//        }
//        double p2p1[3], p4p1[3];
//        for (int i = 0; i < 3; ++i) {
//            p2p1[i] = p2[i] - p1[i];
//            p4p1[i] = p4[i] - p1[i];
//        }
//
//
//        double v1v2v3 = vtkMath::AngleBetweenVectors(p1p2, p3p2);//angle at v2
//        double v2v1v4 = vtkMath::AngleBetweenVectors(p2p1,p4p1);//angle at v1
//
//        if (v1v2v3 > v2v1v4) { 
//            vtkIdType triangle1[3] = { v2, v1, v4 };
//            vtkIdType triangle2[3] = { v3, v2, v4 };
//
//            cells->InsertNextCell(3, triangle1);
//            cells->InsertNextCell(3, triangle2);
//
//
//        }
//        else
//        { 
//            vtkIdType triangle1[3] = { v2, v1, v3 };
//            vtkIdType triangle2[3] = { v1, v4, v3 };
//
//            cells->InsertNextCell(3, triangle1);
//            cells->InsertNextCell(3, triangle2);
//
//        } 
//    }
//     
//    return true;
//}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Create1Quad(vtkPolyData* mesh, vtkIdType e1, vtkIdType e2, double d) {
    vtkPoints* points = mesh->GetPoints();
    if (!points) {
        return false;
    }

    vtkIdType v1 = mesh->GetCell(e1)->GetPointId(0);
    vtkIdType v2 = mesh->GetCell(e1)->GetPointId(1);
    vtkIdType v3 = mesh->GetCell(e2)->GetPointId(0);
    vtkIdType v4 = mesh->GetCell(e2)->GetPointId(1);

    // Calculate vectors for each edge
    double p1[3], p2[3], p3[3], p4[3];

    for (int i = 0; i < 3; ++i) {
        p1[i] = points->GetPoint(v1)[i];
        p2[i] = points->GetPoint(v2)[i];
        p3[i] = points->GetPoint(v3)[i];
        p4[i] = points->GetPoint(v4)[i];
    }

    // Calculate the centroid of the quad
    double centroid[3];
    for (int i = 0; i < 3; ++i) {
        centroid[i] = (p1[i] + p2[i] + p3[i] + p4[i]) / 4.0;
    }

    // Calculate the cross product of two vectors formed by edges
    double crossProduct[3];
    vtkMath::Subtract(p2, p1, crossProduct);
    vtkMath::Subtract(p4, p3, crossProduct);

    // Check the direction of the cross product
    bool clockwise = vtkMath::Dot(centroid, crossProduct) < 0;

    vtkSmartPointer<vtkCellArray> cells = mesh->GetPolys();

    if (clockwise) {
        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
        quad->GetPointIds()->SetId(0, v1);
        quad->GetPointIds()->SetId(1, v2);
        quad->GetPointIds()->SetId(2, v3);
        quad->GetPointIds()->SetId(3, v4);
        cells->InsertNextCell(quad);
    }
    else {
        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
        return false;
        quad->GetPointIds()->SetId(0, v2);
        quad->GetPointIds()->SetId(1, v1);
        quad->GetPointIds()->SetId(2, v4);
        quad->GetPointIds()->SetId(3, v3);
        cells->InsertNextCell(quad);
    }


    return true;
}

//void EnsureOutwardNormals(vtkPolyData* mesh) {
//    // Compute normals
//    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
//    normalGenerator->SetInputData(mesh);
//    normalGenerator->ComputePointNormalsOff();
//    normalGenerator->ComputeCellNormalsOn();
//    normalGenerator->Update();
//
//    // Get the computed normals
//    vtkDataArray* normals = normalGenerator->GetOutput()->GetCellData()->GetNormals();
//
//    // Loop through each triangle
//    vtkSmartPointer<vtkCellArray> cells = mesh->GetPolys();
//    cells->InitTraversal();
//    vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
//    while (cells->GetNextCell(cellPoints)) {
//        double normal[3];
//        normals->GetTuple(cells->GetTraversalLocation(), normal);
//
//        // Check if the normal points inward
//        if (normal[2] < 0) { // Assuming the cylinder is oriented along the z-axis
//            // Swap the vertices of the triangle
//            vtkIdType tmp = cellPoints->GetId(0);
//            cellPoints->SetId(0, cellPoints->GetId(2));
//            cellPoints->SetId(2, tmp);
//            mesh->ReplaceCell(cells->GetTraversalLocation(), cellPoints);
//        }
//    }
//}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

    // deep copy for initialization
    vtkNew<vtkPolyData> mesh;
    mesh->ShallowCopy(input);

    // Create a scalar array to store the slice numbers for points
// Include the necessary header


    //------------------------------------------------------------------------------------------------------------------------------------------------
// Create a scalar array to store the slice numbers for points
    vtkSmartPointer<vtkDoubleArray> v_sliceno = vtkSmartPointer<vtkDoubleArray>::New();
    v_sliceno->SetName("PointSliceNumber");

    // Calculate and assign the slice numbers for points
    vtkIdType numPoints = mesh->GetNumberOfPoints();
    vtkIdType numVertices = mesh->GetNumberOfVerts();

    vtkIdType numEdges = mesh->GetNumberOfLines();
    double prevZ = mesh->GetPoint(0)[2];
    int integral = 1;
    for (vtkIdType i = 0; i < numPoints; ++i) {
        double z = mesh->GetPoint(i)[2];
        if (z > prevZ) {
            integral++;
        }
        v_sliceno->InsertNextValue(integral);
        prevZ = z;
    }
    int slices_begin = 1;
    int slices_end = v_sliceno->GetValue(numPoints - 1);
    //------------------------------------------------------------------------------------------------------------------------------------------------
    //vtkIdType numEdges = 0;
    vtkCellArray* lines = mesh->GetLines();
    logfile << "\n num of points  = " << numPoints;
    logfile << "\n num of Vertices = " << numVertices;
    logfile << "\n num of edges   = " << numEdges;
    //------------------------------------------------------------------------------------------------------------------------------
        // Create a scalar array to store the slice numbers for edges
    vtkSmartPointer<vtkIntArray> e_sliceno = vtkSmartPointer<vtkIntArray>::New();
    e_sliceno->SetName("EdgeSliceNumber");
    //set e_sliceno for all edges //copy from vertices 
    e_sliceno->SetNumberOfValues(3 * numEdges);
    for (vtkIdType e = 0; e < 3 * numEdges; ++e) { e_sliceno->SetValue(e, -1); }
    //------------------------------------------------------------------------------------------------------------------------------------------------ 
    for (vtkIdType e = 0; e < numEdges; e++) { 
                vtkIdType v1= mesh->GetCell(e)->GetPointId(0);
                vtkIdType v2= mesh->GetCell(e)->GetPointId(1);

                double v1_slice = v_sliceno->GetValue(v1);
                double v2_slice = v_sliceno->GetValue(v2);
                if (v1_slice == v2_slice) {
                    e_sliceno->SetValue(e, v1_slice);
                }
                else {
                    e_sliceno->SetValue(e, -1);
                }
                //logfile << "\n dd  e_sliceno[" << e << "] = " << e_sliceno->GetValue(e);
               // logfile << "\t e"<<e<<" = >  " << v1 << "___" << v2;
            
        

    }

    //------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
     /*  for (vtkIdType e = 0; e < numEdges; e++) {
           vtkIdType i = e * 3;
                 logfile << "\n   e_sliceno[" << i << "] = " << e_sliceno->GetValue(i);
    }*/


    // Populate the vectors with edge indices based on their slice numbers
    std::vector<std::vector<vtkIdType>> edges_slice_s(9000);
    std::vector<double> slice_avg_width(9000, 0);//to store average 
    std::vector<int> slice_validpairs_count(9000, 0);//to store average 
    

    // Counter variable to keep track of edge ID within each slice
   // std::vector<size_t> edge_counter(100, 0);

    for (vtkIdType e = 0; e < numEdges; e++) { 
        vtkIdType slice_number = e_sliceno->GetValue(e);

        if (slice_number >= slices_begin && slice_number <= slices_end) {
            edges_slice_s[slice_number - slices_begin].push_back(e);
        }
    }

      //for (vtkIdType s = slices_begin; s <= slices_end; s++) {
      //   // logfile << "Slice " << s << " edges:" << std::endl;

      //    const auto& edges = edges_slice_s[s - slices_begin];
      //    logfile <<"S_"<<s<< "   edges = " << edges.size();
      //    for (size_t i = 0; i < edges.size(); ++i) {
      //      //  logfile << "e_" << edges[i];
      //    }
      //    logfile << std::endl;
      //}
      //------------------------------------------------------------------------------------------------------------------------------------------------------
      //e_up start of forward labelling 
      // Create a vector to store the nearest edge indices from the next slice
    std::vector<vtkIdType> e_up(numEdges, -1); // Initialize with -1 indicating no nearest edge found
    std::vector<vtkIdType> e_btm(numEdges, -1); // Initialize with -1 indicating no nearest edge found
    std::vector<double> e_upd(numEdges, std::numeric_limits<double>::max()); // distance to the nearest edge to avoid conflict with two eges from slice s

    // Iterate over each edge at each slice to find the nearest edge from the next slice
    for (vtkIdType s = slices_begin; s <= slices_end; s++) {
        const auto& edges_btm_slice = edges_slice_s[s - slices_begin];
        const auto& edges_up_slice = edges_slice_s[s - slices_begin + 1];
        double min_dist = std::numeric_limits<double>::max();

        for (vtkIdType e : edges_btm_slice) {
            vtkIdType min_dist_edge = -1; // Initialize with -1 indicating no nearest edge found
            min_dist = std::numeric_limits<double>::max();
            // Iterate over edges in the next slice to find the nearest one
            for (vtkIdType e2 : edges_up_slice) {
                double dist = dist_bw_2es(mesh, e, e2);
                if (dist < min_dist) {
                    min_dist = dist;
                    min_dist_edge = e2;
                }
            }
            e_up[e] = min_dist_edge;
            if (e_upd[min_dist_edge] > min_dist)//valid pair found 
            {
                e_upd[min_dist_edge] = min_dist;
                e_btm[min_dist_edge] = e;
                slice_avg_width[s]+= min_dist;
                slice_validpairs_count[s]++;
            }//else it will be a specail case one upper edge has more than one btm edges in nearest --- so we need one triangle instead of 2 

          // logfile <<"\nSlice_"<<s<< " =  e_up["<<e<<"]" << e_up[e] << "\t min_dist" << min_dist;
        }

        if (slice_validpairs_count[s]>0)
        {
            slice_avg_width[s] /= slice_validpairs_count[s];
        }
    }

    //e_up end --- forward labelling done 
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // now it is time to triangulate each pair of edges ebtm and eup //main loop 
    double zipper_alfa = this->ZipperAlpha;
    vtkNew<vtkCellArray> cells;
     
    for (vtkIdType s = slices_begin; s < slices_end; s++) {
        const auto& edges_btm_slice = edges_slice_s[s - slices_begin];

        for (vtkIdType e1 : edges_btm_slice) {
            // Iterate over edges in the next slice to find the nearest one
            vtkIdType e2 = e_up[e1];

            // Check if e1 and e2 form a valid pair
            if (e2 != -1 && e_btm[e2] == e1) {
                // Create1Quad(output, e1, e2, 2);
                if (e_upd[e2] < zipper_alfa * slice_avg_width[s])
                Create2Triangles(mesh, cells, e1, e2, 1);


            }
            else {
                // Handle the case where e1 and e2 do not form a valid pair
                if (e2 != -1 && e_btm[e2] != -1) {//bottom to up single trinagles filling 
                    //if(1.2*e_upd[e2]>dist_bw_2es(mesh, e1, e2))
                    if (zipper_alfa * slice_avg_width[s] > dist_bw_2es(mesh, e1, e2)) 
                    {
                    Create1Triangle(mesh, cells, e1, e2, e_btm[e2], 1);
                    }

                }

            }
        }


        //--------------------------------------------------------------- top to bottom single triangle filling------------------------------

        const auto& edges_up_slice = edges_slice_s[s - slices_begin + 1];
        for (vtkIdType e2 : edges_up_slice) {
            if (e_btm[e2] == -1)
            {

                vtkIdType min_dist_edge = -1; // Initialize with -1 indicating no nearest edge found
                double min_dist = std::numeric_limits<double>::max();
                // Iterate over edges in the next slice to find the nearest one
                for (vtkIdType e1 : edges_btm_slice) {
                    double dist = dist_bw_2es(mesh, e1, e2);
                    if (dist < min_dist) {
                        min_dist = dist;
                        min_dist_edge = e1;
                    }
                }
                e_btm[e2] = min_dist_edge;
                if (e_upd[e2]>min_dist)
                {
                    e_upd[e2] = min_dist;
                }
                vtkIdType e2adj = e_up[e_btm[e2]];
                if (e2adj != -1) {
                   // if ( e_upd[e2] < 1.2 * dist_bw_2es(mesh, e2adj, e_btm[e2adj]))

                        if (e_upd[e2] <zipper_alfa * slice_avg_width[s])
                        {
                            Create1Triangle(mesh, cells, e2, e_btm[e2], e2adj, 1);
                        }
                }


            }

        }

        //---------------------------------------------------------------repppp

        logfile << "\nSlice_" << s;
    }
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
        // Close the text file
    logfile.close();
     
    double colocate_epsilon = 0.1; 
    vtkNew<vtkPolyData> final;
    vtkNew<vtkPoints> points;
    for (auto k = 0; k < mesh->GetNumberOfPoints(); k++)
        points->InsertNextPoint(mesh->GetPoint(k));

    final->SetPoints(points);
    final->SetPolys(cells);
     

   // repair_colocate_vertices(final, colocate_epsilon);

     
    vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    output->ShallowCopy(final);


    return 1;
}

int vtkZipperTriangulation::FillInputPortInformation(int port, vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}
