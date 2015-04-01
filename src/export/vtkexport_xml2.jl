using LightXML
using Codecs

abstract AbstractVTKXML

#=
type VTKXMLCompressed <: AbstractVTKXML
    buffer::IOBuffer
    compressed
end

push!(vtkxml::AbstractVTKXML, data) = write(vtkxml, buffer, data)

function get_header(vtkxml::VTKXMLCompressed)
    return UInt32[1, buffer.size, buffer.size, length(vtkxml.compressed)]
end

get_header(vtkxml::VTKXMLUncompressed) = UInt32(buffer.size)


function getsize(compressed::Bool) =

function get_utf8_size()

function get_
    bytestring(encode(Base64, reinterpret(UInt8, size_block)))))

function get_
=#

typealias Vtkd Dict{ASCIIString, ASCIIString}

function write_VTKXML_section{T <: FESection}(filename, nodes::Vector{FENode2},
                                              section::T, binary::Bool=true, compress::Bool=true)
    if !binary && compress
        error("Can only compress when using Binary format")
    end

    if binary
        const VTK_FORMAT = "binary"
    else
        const VTK_FORMAT = "ascii"

    iobuf = IOBuffer()
    xdoc = XMLDocument()
    xroot = create_root(xdoc, "VTKFile")
    set_attribute(xroot, "type", "UnstructuredGrid")
    set_attribute(xroot, "version", "0.1")
    set_attribute(xroot, "byte_order", "LittleEndian")

    if compress
        set_attribute(xroot, "compressor", "vtkZLibDataCompressor")
    end


    xgrid = new_child(xroot, "UnstructuredGrid")

    xpiece = new_child(xgrid, "Piece")
    set_attribute(xpiece, "NumberOfPoints", fp.nodes)


    ncells = length(section.elements)
    set_attribute(xpiece, "NumberOfCells", ncells)

    # Points
    xpoints = new_child(xpiece, "Points")

    # Coordinates for points
    xcoords = new_child(xpoints, "DataArray")
    set_attribute(xcoords, Vtkd("type" => "Float64", "name" => "Points",
                                "format" => VTK_FORMAT), "NumberOfComponents" => "3")

    # Write points binary
    if binary
        # If we are not compressing we should first write the
        # Base64 encoded UInt32 value of the number of bytes
        if !compress
            nodecoords_size = reinterpret(UInt8, [UInt32(3*length(nodes))])
            add_text(xcoords, bytestring(encode(Base64, nodecoords_size)))
        end

        for node in nodes
            coords = get_coord(node)
            for c in coords
                write(iobuf, c)
            end
        end

        if !compress
            nodecoords_b64 = encode(Base64, takebuf_array(iobuf))
            add_text(bytestring(nodescoords_b64))

        else
            uncomp_size = iobuf.size
            encoded_data = encode(Zlib, takebuf_array(iobuf))
            compressed_size = length(encoded_data)
            size_block = map(UInt32, [1; uncomp_size; uncomp_size; compressed_size])
            add_text(xcoords, bytestring(encode(Base64, reinterpret(UInt8, size_block)))))
            add_text(xcoords, bytestring(encode(Base64, encoded_data)))
        end

    else # Write points ASCII
        for node in nodes
            coords = get_coord(node)
            for c in coords
                print(iobuf, c, " ")
            end
        end
        add_text(xcoords, takebuf_string(iobuf))
    end

    #### Cells ####
    xcells = new_child(xpiece, "Cells")

    # Cell locations
    xcellconn = new_child(xcells, "DataArray")
    set_attribute(xcellconn, "type", "Int64")
    set_attribute(xcellconn, "Name", "connectivity")
    set_attribute(xcellconn, "format", VTK_FORMAT)

    # Cell location data
    write(iobuf, "\n")
    write(iobuf, "0 1 2 \n")


    xcell_offsets = new_child(xcells, "DataArray")
    set_attribute(xcell_offsets, "type", "Int64")
    set_attribute(xcell_offsets, "Name", "offsets")
    set_attribute(xcell_offsets, "format", VTK_FORMAT)
    nverts = length(ele_1.vertices)
    offsets = collect(nverts:nverts:nverts*length(section.elements))



    xcell_types = new_child(xcells, "DataArray")
    set_attribute(xcell_types, "type", "UInt8")
    set_attribute(xcell_types, "Name", "types")
    set_attribute(xcell_types, "format", VTK_FORMAT)
    ele_1 = section.elements[1]
    cell_types = [get_vtk_num(get_geotype(ele1)) for i in 1:length(section.elements)]


    #### Cell data (dummy) ####
    xcell_data = new_child(xpiece, "CellData")
    xcol = new_child(xcell_data, "DataArray")
    set_attribute(xcol, "Name", "Stress")
    set_attribute(xcol, "type", "Float64")
    set_attribute(xcol, "format", VTK_FORMAT)



    for element in section.elements

        add_text



    write(iobuf, "\n")
    write(iobuf, "5 \n")
    add_text(xcell_types, takebuf_string(iobuf))


    #### Data at Points ####
    xpoint_data = new_child(xpiece, "PointData")

    # Points
    xdisp = new_child(xpoint_data, "DataArray")
    set_attribute(xdisp, "Name", "Displacement")
    set_attribute(xdisp, "NumberOfComponents", "3")
    set_attribute(xdisp, "type", "Float64")

    #=
    if binary
        set_attribute(xdisp, "format", "binary")
        data = vec([0.0 0.0 0.0 0.0 0.0 0.0 10.0 10.0 10.0])
        #write(iobuf, bytestring(encode(Base64, reinterpret(UInt8, [UInt32(8*9)]))))
        write(iobuf, UInt32(8*9))
        write(iobuf, data)
        #write(iobuf, encode(Base64, bytestring(reinterpret(UInt8, data))))
        #add_text(xdisp, takebuf_string(iobuf))
        add_text(xdisp, bytestring(encode(Base64, takebuf_string(iobuf))))
    else
        =#
        set_attribute(xdisp, "format", "ascii")
        write(iobuf, "\n")
        write(iobuf, "0 0 0 0 0 0 10 10 10\n")
        add_text(xdisp, takebuf_string(iobuf))
    #end

    write(iobuf, "\n")
    write(iobuf, " 1\n")
    add_text(xcol, takebuf_string(iobuf))

    save_file(xdoc, filename)
end
