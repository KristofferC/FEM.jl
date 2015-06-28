module VTKExportMod

using FEM

using Zlib
using LightXML
using Codecs

import Base.push!
import FEM: write_data, get_coord, get_geotype, get_ncomponents, get_cell_data,
            get_displacement, get_field, get_galphas, DofVals
export VTKExporter, set_binary!, set_compress!



type VTKExporter <: AbstractDataExporter
    fields::Vector{DataType}
    binary::Bool
    compress::Bool
end
fields(vtkexp::VTKExporter) = vtkexp.fields
VTKExporter() = VTKExporter(Vector{DataType}[], true, false)

function set_binary!(vtkexp::VTKExporter, val::Bool)
    vtkexp.binary = val
    if !val
        vtkexp.compress = false
    end
end

function set_compress!(vtkexp::VTKExporter, val::Bool)
    if !vtkexp.binary && val
        warning("Can only compress when VTKEXporter is in binary mode")
        return
    else
        vtkexp.compress = val
    end
end
push!(vtkexp::VTKExporter, field::DataType) = push!(vtkexp.fields, field)


#Long ass names in spirit of the VTK library :)
abstract AbstractVTKXMLWriter

abstract AbstractVTKXMLBinaryWriter <: AbstractVTKXMLWriter

function add_data!(vtkw::AbstractVTKXMLBinaryWriter, data)
    write(vtkw.buffer, data)
end

type VTKXMLBinaryCompressedWriter <: AbstractVTKXMLBinaryWriter
    buffer::IOBuffer
end
VTKXMLBinaryCompressedWriter() = VTKXMLBinaryCompressedWriter(IOBuffer())

function write_data!(vtkw::VTKXMLBinaryCompressedWriter, xmlele::XMLElement)
    uncompressed_size = vtkw.buffer.size
    buff = takebuf_array(vtkw.buffer)
    #compressed_data = encode(Zlib, buff, 5)
    compressed_data = compress(buff, 5)
    compressed_size = length(compressed_data)
    header = UInt32[1, uncompressed_size, uncompressed_size, compressed_size]
    header_binary = bytestring(encode(Base64, reinterpret(UInt8, header)))
    data_binary = bytestring(encode(Base64, compressed_data))
    add_text(xmlele, header_binary)
    add_text(xmlele, data_binary)
end

type VTKXMLBinaryUncompressedWriter <: AbstractVTKXMLBinaryWriter
    buffer::IOBuffer
end
VTKXMLBinaryUncompressedWriter() = VTKXMLBinaryUncompressedWriter(IOBuffer())

function write_data!(vtkw::VTKXMLBinaryUncompressedWriter, xmlele::XMLElement)
    uncompressed_size = vtkw.buffer.size
    uncompressed_data = takebuf_array(vtkw.buffer)
    header = UInt32[uncompressed_size]
    header_binary = bytestring(encode(Base64, reinterpret(UInt8, header)))
    data_binary = bytestring(encode(Base64, uncompressed_data))
    add_text(xmlele, header_binary)
    add_text(xmlele, data_binary)
end


type VTKXMLASCIIWriter <: AbstractVTKXMLWriter
    buffer::IOBuffer
end
VTKXMLASCIIWriter() = VTKXMLASCIIWriter(IOBuffer())


function add_data!(vtkw::VTKXMLASCIIWriter, data)
    print(vtkw.buffer, data, " ")
end

function add_data!{T <: AbstractArray}(vtkw::VTKXMLASCIIWriter, data::T)
    for comp in data
        print(vtkw.buffer, comp, " ")
    end
end

function write_data!(vtkw::VTKXMLASCIIWriter, xmlele::XMLElement)
    add_text(xmlele, takebuf_string(vtkw.buffer))
end

get_vtk_num(::GeoTrig) = 5
get_vtk_num(::Type{GeoTrig}) = 5
get_vtk_num(::GeoQTrig) = 22
get_vtk_num(::Type{GeoQTrig}) = 22
get_vtk_num(::GeoQuad) = 9
get_vtk_num(::Type{GeoQuad}) = 9
get_vtk_num(::GeoTetra) = 10
get_vtk_num(::Type{GeoTetra}) = 10


function write_data(fp::FEProblem, vtkexp::VTKExporter, tstep::Int)

    if vtkexp.binary
        if vtkexp.compress
            writer = VTKXMLBinaryCompressedWriter()
        else
            writer = VTKXMLBinaryUncompressedWriter()
        end
    else
        writer = VTKXMLASCIIWriter()
    end

    for section in fp.sections
        _write_VTKXML_section(fp.name, fp.nodes, fp.dof_vals, section, vtkexp, writer, tstep)
    end
end


function _write_VTKXML_section(filename::ASCIIString, nodes::Vector{FENode2}, dof_vals::DofVals,
                               section::FESection, vtkexp:: VTKExporter, vtkw::AbstractVTKXMLWriter, tstep::Int)

    if vtkexp.binary
        const VTK_FORMAT = "binary"
    else
        const VTK_FORMAT = "ascii"
    end

    xdoc = XMLDocument()
    xroot = create_root(xdoc, "VTKFile")
    set_attribute(xroot, "type", "UnstructuredGrid")
    set_attribute(xroot, "version", "0.1")
    set_attribute(xroot, "byte_order", "LittleEndian")

    if vtkexp.compress
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
    set_attribute(xcoords, "type", "Float64")
    set_attribute(xcoords, "name", "Points")
    set_attribute(xcoords, "format", VTK_FORMAT)
    set_attribute(xcoords, "NumberOfComponents", "3")

    for node in nodes
        coords = get_coord(node)
        for c in coords
            add_data!(vtkw, c)
        end
    end
    write_data!(vtkw, xcoords)

    # Cells
    xcells = new_child(xpiece, "Cells")

    # Cell connectivity
    xcellconn = new_child(xcells, "DataArray")
    set_attribute(xcellconn, "type", "Int64")
    set_attribute(xcellconn, "Name", "connectivity")
    set_attribute(xcellconn, "format", VTK_FORMAT)

    for element in section.elements
        for vertex in element.vertices
            add_data!(vtkw, vertex-1)
        end
    end
    write_data!(vtkw, xcellconn)


    # Cell location data
    xcell_offsets = new_child(xcells, "DataArray")
    set_attribute(xcell_offsets, "type", "Int64")
    set_attribute(xcell_offsets, "Name", "offsets")
    set_attribute(xcell_offsets, "format", VTK_FORMAT)
    nverts = length(section.elements[1].vertices)
    offsets = collect(nverts:nverts:nverts*length(section.elements))
    add_data!(vtkw, offsets)
    write_data!(vtkw, xcell_offsets)


    xcell_types = new_child(xcells, "DataArray")
    set_attribute(xcell_types, "type", "UInt8")
    set_attribute(xcell_types, "Name", "types")
    set_attribute(xcell_types, "format", VTK_FORMAT)
    ele_vtknum = get_vtk_num(get_geotype(section.elements[1]))
    cell_types = UInt8[ele_vtknum for _ in 1:length(section.elements)]
    add_data!(vtkw, cell_types)
    write_data!(vtkw, xcell_types)


    # Cell data
    xcell_data = new_child(xpiece, "CellData")
    for field in fields(vtkexp)
        # Separated this into its own function for type stability w.r.t. fields
        xcellfield = new_child(xcell_data, "DataArray")
        write_celldata_field!(xcellfield, field, section, vtkw, VTK_FORMAT)
    end


    # Point data
    xpoint_data = new_child(xpiece, "PointData")

    # Points
    xdisp = new_child(xpoint_data, "DataArray")
    set_attribute(xdisp, "Name", "Displacement")
    set_attribute(xdisp, "NumberOfComponents", "3")
    set_attribute(xdisp, "type", "Float64")
    set_attribute(xdisp, "format", VTK_FORMAT)
    for node in nodes
        disp = get_displacement(node, dof_vals)
        for i in disp
            add_data!(vtkw, i)
        end
    end
    write_data!(vtkw, xdisp)




    println("Finished writing data")

    save_file(xdoc, string(filename, "_$(tstep).vtu"))

    #free(xdoc) # Does this free everything?
end

function write_celldata_field!{T <: AbstractField}(xcellfield::XMLElement, field::Type{T},
                               section::FESection, vtkw::AbstractVTKXMLWriter, VTK_FORMAT::ASCIIString)
    set_attribute(xcellfield, "NumberOfComponents", get_ncomponents(field))
    set_attribute(xcellfield, "Name", string(field))
    set_attribute(xcellfield, "type", "Float64")
    set_attribute(xcellfield, "format", VTK_FORMAT)
    for elem in section.elements
        add_data!(vtkw, get_cell_data(elem, field))
    end
    write_data!(vtkw, xcellfield)
end

end # module
