

# Export mesh to a VTK 1.0 file as an unstructured grid.
function vtkexportmesh (theFile::String,Connectivity,Points, Cell_type;
                        vectors=nothing,vectors_name ="vectors",
                        scalars=nothing, scalars_name ="scalars",binary = false)



    fid=open(theFile,"w");
    if (fid==-1)
        error (["Could not open " * theFile])
        return nothing
    end

    print(fid,"# vtk DataFile Version 1.0\n");
    if (!binary)
        print(fid,"ASCII\n");
    else
        print(fid,"BINARY\n");
    end
    print(fid,"\n");
    print(fid,"DATASET UNSTRUCTURED_GRID\n");
    print(fid,"POINTS ", size(X,1), " double\n");
    if (!binary)
        for node in mesh.nodes
                print(fid,node.coords[1],node.coords[2], node.coords[3],"\n")
        end
    else
        fwrite(fid,cast(X,"double"),"double","n");
    end
    print(fid,"\n");
    print(fid,"\n");

    n_els = length(mesh.elements)
    print(fid,"CELLS ",n_els," ",(n_els*(size(Connectivity,2)+1)),"\n");
    if (!binary)
        for i= 1:size(Connectivity, 1)
            print(fid,size(Connectivity,2)," ");
            for j= 1:size(Connectivity,2)-1
                print(fid,Connectivity[i,j]-1," ");
            end
            print(fid,Connectivity[i,end]-1,"\n");
        end
    else
        fwrite(fid,cast([zeros(size(f,1),1)+size(f,2),f-1],"int32"),"int32","n");
    end
    print(fid,"\n");
    print(fid,"\n");
    print(fid,"CELL_TYPES ",size(Connectivity,1),"\n");
    if (!binary)
        for i= 1:size(Connectivity,1)
            print(fid,Cell_type,"\n");
        end
    else
        fwrite(fid,cast(ctype,"int32"),"int32","n");
    end

    did_point_data=false
    if (scalars!=nothing) && (length(scalars)==size(Points,1))
        did_point_data=true
        print(fid,"POINT_DATA ",length(scalars),"\n");
        print(fid,"SCALARS ",scalars_name," double\n");
        print(fid,"LOOKUP_TABLE default\n");
        if (!binary)
            for j= 1:length(scalars)
                print(fid,scalars[j],"\n");
            end
        else
            fwrite(fid,cast(scalars,"double"),"double","n");
        end
    end

    if vectors!=nothing
        if (!did_point_data)
            print(fid,"POINT_DATA ",size(vectors,1),"\n");
        end
        did_point_data=true
        print(fid,"VECTORS ",vectors_name," double\n");
        #print(fid,"LOOKUP_TABLE default\n");
        X=vectors
        if size(vectors, 2)<3
            X=   zeros(size (vectors,1),3)
        end
        X[:,1:size(vectors,2)] = vectors
        if (!binary)
            for j= 1:size(X,1)
                k=1;
                print(fid,X[j,k]);
                for k=2:size(X,2)
                    print(fid," ",X[j,k]);
                end
                print(fid,"\n");
            end
        else
            fwrite(fid,cast(X,"double"),"double","n");
        end
    end

    if (scalars!=nothing) && (length(scalars)==size(Connectivity,1))
        did_point_data=true
        print(fid,"CELL_DATA ",length(scalars),"\n");
        print(fid,"SCALARS ",scalars_name," double\n");
        print(fid,"LOOKUP_TABLE default\n");
        if (!binary)
            for j= 1:length(scalars)
                print(fid,scalars[j],"\n");
            end
        else
            fwrite(fid,cast(scalars,"double"),"double","n");
        end
    end

    print(fid,"\n");
    fid=close(fid);
    return true
end
export vtkexportmesh


