using LightXML
using Codecs

abstract AbstractVTKXML

abstract AbstractVTKXMLBinary <: AbstractVTKXML

function add_data!{T <: AbstractVTKXMLBinary}(vtkxml::T, data)
    write(vtkxml.buffer, data)
end

type VTKXMLBinaryCompressed <: AbstractVTKXMLBinary
    buffer::IOBuffer
end
VTKXMLBinaryCompressed() = VTKXMLBinaryCompressed(IOBuffer())

function write_data!(vtkxml::VTKXMLBinaryCompressed, xmlele::XMLElement)
    uncompressed_size = vtkxml.buffer.size
    compressed_data = encode(Zlib, takebuf_array(vtkxml.buffer))
    compressed_size = length(compressed_data)
    header = UInt32[1, uncompressed_size, uncompressed_size, compressed_size]
    header_binary = bytestring(encode(Base64, reinterpret(UInt8, header)))
    data_binary = bytestring(encode(Base64, compressed_data))
    add_text(xmlele, header_binary)
    add_text(xmlele, data_binary)
end

type VTKXMLBinaryUncompressed <: AbstractVTKXMLBinary
    buffer::IOBuffer
end
VTKXMLBinaryUncompressed() = VTKXMLBinaryUncompressed(IOBuffer())

function write_data!(vtkxml::VTKXMLBinaryUncompressed, xmlele::XMLElement)
    uncompressed_size = vtkxml.buffer.size
    uncompressed_data = takebuf_array(vtkxml.buffer)
    header = UInt32[uncompressed_size]
    header_binary = bytestring(encode(Base64, reinterpret(UInt8, header)))
    data_binary = bytestring(encode(Base64, uncompressed_data))
    add_text(xmlele, header_binary)
    add_text(xmlele, data_binary)
end


type VTKXMLASCII <: AbstractVTKXML
    buffer::IOBuffer
end
VTKXMLASCII() = VTKXMLASCII(IOBuffer())


function add_data!(vtkxml::VTKXMLASCII, data)
    print(vtkxml.buffer, data, " ")
end

function add_data!{T <: AbstractArray}(vtkxml::VTKXMLASCII, data::T)
    for comp in data
        print(vtkxml.buffer, comp, " ")
    end
end

function write_data!(vtkxml::VTKXMLASCII, xmlele::XMLElement)
    add_text(xmlele, takebuf_string(vtkxml.buffer))
end

function write_VTKXML(filename::ASCIIString, fp::FEProblem,
                      binary::Bool=true, compress::Bool=true)
    for section in fp.sections
        write_VTKXML_section(filename, fp.nodes, section, binary, compress)
    end
end

function write_VTKXML_section{T <: FESection}(filename::ASCIIString, nodes::Vector{FENode2},
                             section::T, binary::Bool=true, compress::Bool=false)
    if binary
        if compress
            _write_VTKXML_section(filename, nodes, section, binary, compress, VTKXMLBinaryCompressed())
        else
            _write_VTKXML_section(filename, nodes, section, binary, compress, VTKXMLBinaryUncompressed())
        end
    else
        _write_VTKXML_section(filename, nodes, section, binary, compress, VTKXMLASCII())
    end
end

typealias Vtkd Dict{ASCIIString, ASCIIString}

function _write_VTKXML_section{T <: FESection, P <: AbstractVTKXML}(filename::ASCIIString, nodes::Vector{FENode2},
                                              section::T, binary::Bool, compress::Bool, vtkxml::P)
    if !binary && compress
        error("Can only compress when using Binary format")
    end

    if binary
        const VTK_FORMAT = "binary"
    else
        const VTK_FORMAT = "ascii"
    end

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
    set_attribute(xpiece, "NumberOfPoints", length(nodes))


    ncells = length(section.elements)
    set_attribute(xpiece, "NumberOfCells", length(section.elements))

    # Points
    xpoints = new_child(xpiece, "Points")

    # Coordinates for points
    xcoords = new_child(xpoints, "DataArray")
    set_attributes(xcoords, Vtkd("type" => "Float64", "name" => "Points",
                                "format" => VTK_FORMAT, "NumberOfComponents" => "3"))

    for node in nodes
        coords = get_coord(node)
        for c in coords
            add_data!(vtkxml, c)
        end
    end

    write_data!(vtkxml, xcoords)

    #### Cells ####

    xcells = new_child(xpiece, "Cells")


    # Cell locations
    xcellconn = new_child(xcells, "DataArray")
    set_attribute(xcellconn, "type", "Int64")
    set_attribute(xcellconn, "Name", "connectivity")
    set_attribute(xcellconn, "format", VTK_FORMAT)

    for element in section.elements
        for vertex in element.vertices
            add_data!(vtkxml, vertex-1)
        end
    end
    write_data!(vtkxml, xcellconn)


    # Cell location data


    xcell_offsets = new_child(xcells, "DataArray")
    set_attribute(xcell_offsets, "type", "Int64")
    set_attribute(xcell_offsets, "Name", "offsets")
    set_attribute(xcell_offsets, "format", VTK_FORMAT)
    nverts = length(section.elements[1].vertices)
    offsets = collect(nverts:nverts:nverts*length(section.elements))
    add_data!(vtkxml, offsets)
    write_data!(vtkxml, xcell_offsets)


    xcell_types = new_child(xcells, "DataArray")
    set_attribute(xcell_types, "type", "UInt8")
    set_attribute(xcell_types, "Name", "types")
    set_attribute(xcell_types, "format", VTK_FORMAT)
    ele_vtknum = get_vtk_num(get_geotype(section.elements[1]))
    cell_types = UInt8[ele_vtknum for i in 1:length(section.elements)]
    add_data!(vtkxml, cell_types)
    write_data!(vtkxml, xcell_types)

    #### Cell data (dummy) ####
    #xcell_data = new_child(xpiece, "CellData")
    #xcol = new_child(xcell_data, "DataArray")
    #set_attribute(xcol, "Name", "Stress")
    #set_attribute(xcol, "type", "Float64")
    #set_attribute(xcol, "format", VTK_FORMAT)


 #=

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

        set_attribute(xdisp, "format", "ascii")
        write(iobuf, "\n")
        write(iobuf, "0 0 0 0 0 0 10 10 10\n")
        add_text(xdisp, takebuf_string(iobuf))
    #end

    write(iobuf, "\n")
    write(iobuf, " 1\n")
    add_text(xcol, takebuf_string(iobuf))
    =#

    save_file(xdoc, filename)
end
