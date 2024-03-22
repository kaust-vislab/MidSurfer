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
//...........................................................................................................................  
bool Create2Triangles(vtkPolyData* mesh, vtkIdType e1, vtkIdType e2, double d) {
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

     

    //// Calculate the direction vectors for each edge
    //double dir1[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
    //double dir2[3] = { p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2] }; 
    // Check the direction of edges e1 and e2
    double dotProduct = vtkMath::Dot(dir1, dir2);
    vtkSmartPointer<vtkCellArray> cells = mesh->GetPolys();

    if (dotProduct > 0) {
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
        if (v1v2v4 > v2v1v3) {
            return false;
            vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
            triangle1->GetPointIds()->SetId(0, v1);
            triangle1->GetPointIds()->SetId(1, v2);
            triangle1->GetPointIds()->SetId(2, v3);

            vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
            triangle2->GetPointIds()->SetId(0, v2);
            triangle2->GetPointIds()->SetId(1, v3);
            triangle2->GetPointIds()->SetId(2, v4);

            cells->InsertNextCell(triangle1);
            cells->InsertNextCell(triangle2);

        }
        else {
            vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
            triangle1->GetPointIds()->SetId(0, v1);
            triangle1->GetPointIds()->SetId(1, v2);
            triangle1->GetPointIds()->SetId(2, v4);

            vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
            triangle2->GetPointIds()->SetId(0, v1);
            triangle2->GetPointIds()->SetId(1, v3);
            triangle2->GetPointIds()->SetId(2, v4);

            cells->InsertNextCell(triangle1);
            cells->InsertNextCell(triangle2); 
        }
    }
    else {// opposit direction so dignoal may either be v2v4 OR v1v3

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

        if (v1v2v3 > v2v1v4) {
           // M.facets.create_triangle(v1, v2, v4);
            //M.facets.create_triangle(v2, v4, v3);
            vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
            triangle1->GetPointIds()->SetId(0, v1);
            triangle1->GetPointIds()->SetId(1, v2);
            triangle1->GetPointIds()->SetId(2, v4);

            vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
            triangle2->GetPointIds()->SetId(0, v2);
            triangle2->GetPointIds()->SetId(1, v4);
            triangle2->GetPointIds()->SetId(2, v3);

            cells->InsertNextCell(triangle1);
            cells->InsertNextCell(triangle2);

        }
        else
        {
           // M.facets.create_triangle(v1, v2, v3);
           // M.facets.create_triangle(v1, v3, v4);

            vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
            triangle1->GetPointIds()->SetId(0, v1);
            triangle1->GetPointIds()->SetId(1, v2);
            triangle1->GetPointIds()->SetId(2, v3);

            vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
            triangle2->GetPointIds()->SetId(0, v1);
            triangle2->GetPointIds()->SetId(1, v3);
            triangle2->GetPointIds()->SetId(2, v4);

            cells->InsertNextCell(triangle1);
            cells->InsertNextCell(triangle2);
        } 
    }
     
    return true;
}
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
    int slices_begin = 1;
    int slices_end = v_sliceno->GetValue(numPoints - 1);
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
    e_sliceno->SetNumberOfValues(3 * numEdges);
    for (vtkIdType e = 0; e < 3 * numEdges; ++e) { e_sliceno->SetValue(e, -1); }
    //------------------------------------------------------------------------------------------------------------------------------------------------ 
    for (vtkIdType e = 0; e < numEdges; e++) { 
                vtkIdType v1= output->GetCell(e)->GetPointId(0);
                vtkIdType v2= output->GetCell(e)->GetPointId(1);

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
    std::vector<std::vector<vtkIdType>> edges_slice_s(100);

    // Counter variable to keep track of edge ID within each slice
   // std::vector<size_t> edge_counter(100, 0);

    for (vtkIdType e = 0; e < numEdges; e++) { 
        vtkIdType slice_number = e_sliceno->GetValue(e);

        if (slice_number >= slices_begin && slice_number <= slices_end) {
            edges_slice_s[slice_number - slices_begin].push_back(e);
        }
    }

    /*  for (vtkIdType s = slices_begin; s < slices_end; s++) {
          logfile << "Slice " << s << " edges:" << std::endl;

          const auto& edges = edges_slice_s[s - slices_begin];
          logfile << "Total = " << edges.size();
          for (size_t i = 0; i < edges.size(); ++i) {
              logfile << "e_" << edges[i];
          }
          logfile << std::endl;
      }*/
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
                double dist = dist_bw_2es(output, e, e2);
                if (dist < min_dist) {
                    min_dist = dist;
                    min_dist_edge = e2;
                }
            }
            e_up[e] = min_dist_edge;
            if (e_upd[min_dist_edge] > min_dist)
            {
                e_upd[min_dist_edge] = min_dist;
                e_btm[min_dist_edge] = e;
            }//else it will be a specail case one upper edge has more than one btm edges in nearest --- so we need one triangle instead of 2 

           logfile <<"\nSlice_"<<s<< " =  e_up["<<e<<"]" << e_up[e] << "\t min_dist" << min_dist;
        }
    }

    //e_up end --- forward labelling done 
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // now it is time to triangulate each pair of edges ebtm and eup 

    for (vtkIdType s = slices_begin; s < slices_end; s++) {
        const auto& edges_btm_slice = edges_slice_s[s - slices_begin]; 

        for (vtkIdType e1 : edges_btm_slice) {
            // Iterate over edges in the next slice to find the nearest one
            vtkIdType e2 = e_up[e1];

            // Check if e1 and e2 form a valid pair
            if (e2 != -1 && e_btm[e2] == e1) {
                // Get the vertices of edges e1 and e2 happy d
               // Create1Quad(output, e1, e2, 2);
                Create2Triangles(output, e1, e2, 2);
               // vtkIdType v0 = output->GetCell(e1)->GetPointId(0);
               // vtkIdType v1 = output->GetCell(e1)->GetPointId(1);

               // vtkIdType v2 = output->GetCell(e2)->GetPointId(0);
               // vtkIdType v3 = output->GetCell(e2)->GetPointId(1);

               // // Create a quad from the vertices of e1 and e2  

               // // Now, you can use quad to create a vtkCell or perform further operations 
               //// Create two triangles
               // vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
               // triangle1->GetPointIds()->SetId(0, v0);
               // triangle1->GetPointIds()->SetId(1, v1);
               // triangle1->GetPointIds()->SetId(2, v2);

               // vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
               // triangle2->GetPointIds()->SetId(0, v0);
               // triangle2->GetPointIds()->SetId(1, v2);
               // triangle2->GetPointIds()->SetId(2, v3);

               // // Add triangles to your vtkCellArray// Get the cell array from the output polydata
               // vtkSmartPointer<vtkCellArray> cells = output->GetPolys();

               // cells->InsertNextCell(triangle1);
               // cells->InsertNextCell(triangle2);

            }
            else {
                // Handle the case where e1 and e2 do not form a valid pair
            }
        }
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






    //---------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------
//    The following loop is used to label each edge with e_btm and e_up being its two neares edges from the neibhouring slices  
   //for (vtkIdType s = slices_begin; s < slices_end; s++) {
   //    int count = 0;
   //    double davg = 0.0;  
   //    for (vtkIdType e = 0; e < numEdges; e++) {
   //        vtkIdType i = e * 3;
   //        if (e_sliceno->GetValue(i) != s)
   //            continue;
   //        double mind2 = std::numeric_limits<double>::max();

   //        for (vtkIdType e2 = 0; e2 < numEdges; e2++) {
   //            vtkIdType i2 = e2 * 3;
   //            if (e_sliceno->GetValue(i2) != s+1)
   //                continue;
   //            double de2 = dist_bw_2es(output, i, i2);

   //           // logfile << "\n   Slice # [" << s << "] de2= " << de2;
   //            if (mind2 > de2)
   //            {
   //                mind2 = de2; 
   //               // e_up[e] = e2; 
   //            }

   //        }

   //        davg += mind2;
   //        count++;
   //    }
   //    if(count>0)
   //    davg /= count;

   //    logfile << "\n   Slice # [" << s << "] davg= " << davg;
   //    logfile << "\n   Slice # [" << s << "] count= " << count;

   //}
       //-----------------------------------------------------------------------------------------------------------------------------------------------------------------


        // Close the text file
    logfile.close();




    return 1;
}

int vtkZipperTriangulation::FillInputPortInformation(int port, vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}
