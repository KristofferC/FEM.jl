using LightXML


        doc = XMLDocument()
        xroot = create_root(doc, "VTKFile")
        set_attribute(xroot, "type", "UnstructuredGrid")
        set_attribute(xroot, "version", "0.1")
        set_attribute(xroot, "byte_order", "LittleEndian")

        # Unstructured grid element
        xunstructured_grid = new_child(xroot, "UnstructuredGrid")

        # Piece 0 (only one)
        xpiece = new_child(xunstructured_grid, "Piece")
        set_attribute(xpiece, "NumberOfPoints", 5)
        set_attribute(xpiece, "NumberOfCells", "0")

        ### Points ####
        xpoints = new_child(xpiece, "Points")

        # Point location data
        xpoint_coords = new_child(xpoints, "DataArray")
        set_attribute(xpoint_coords, "type", "Float64")
        set_attribute(xpoint_coords, "format", "ascii")
        set_attribute(xpoint_coords, "NumberOfComponents", "3")
        points.appendChild(point_coords)

        string = self.coords_to_string(x, y, z)
        point_coords_data = doc.createTextNode(string)
        point_coords.appendChild(point_coords_data)

        #### Cells ####
        cells = doc.createElementNS("VTK", "Cells")
        piece.appendChild(cells)

        # Cell locations
        cell_connectivity = doc.createElementNS("VTK", "DataArray")
        cell_connectivity.set_attribute("type", "Int32")
        cell_connectivity.set_attribute("Name", "connectivity")
        cell_connectivity.set_attribute("format", "ascii")
        cells.appendChild(cell_connectivity)

        # Cell location data
        connectivity = doc.createTextNode("0")
        cell_connectivity.appendChild(connectivity)

        cell_offsets = doc.createElementNS("VTK", "DataArray")
        cell_offsets.set_attribute("type", "Int32")
        cell_offsets.set_attribute("Name", "offsets")
        cell_offsets.set_attribute("format", "ascii")
        cells.appendChild(cell_offsets)
        offsets = doc.createTextNode("0")
        cell_offsets.appendChild(offsets)

        cell_types = doc.createElementNS("VTK", "DataArray")
        cell_types.set_attribute("type", "UInt8")
        cell_types.set_attribute("Name", "types")
        cell_types.set_attribute("format", "ascii")
        cells.appendChild(cell_types)
        types = doc.createTextNode("1")
        cell_types.appendChild(types)

        #### Data at Points ####
        point_data = doc.createElementNS("VTK", "PointData")
        piece.appendChild(point_data)

        # Points
        point_coords_2 = doc.createElementNS("VTK", "DataArray")
        point_coords_2.set_attribute("Name", "Points")
        point_coords_2.set_attribute("NumberOfComponents", "3")
        point_coords_2.set_attribute("type", "Float32")
        point_coords_2.set_attribute("format", "ascii")
        point_data.appendChild(point_coords_2)

        string = self.coords_to_string(x, y, z)
        point_coords_2_Data = doc.createTextNode(string)
        point_coords_2.appendChild(point_coords_2_Data)

        # Particle jump vectors
        if len(x_jump) > 0:
            jumps = doc.createElementNS("VTK", "DataArray")
            jumps.set_attribute("Name", "jumps")
            jumps.set_attribute("NumberOfComponents", "3")
            jumps.set_attribute("type", "Float32")
            jumps.set_attribute("format", "ascii")
            point_data.appendChild(jumps)

            string = self.coords_to_string(x_jump, y_jump, z_jump)
            jumpData = doc.createTextNode(string)
            jumps.appendChild(jumpData)

        # Force vectors
        if len(x_force) > 0:
            forces = doc.createElementNS("VTK", "DataArray")
            forces.set_attribute("Name", "forces")
            forces.set_attribute("NumberOfComponents", "3")
            forces.set_attribute("type", "Float32")
            forces.set_attribute("format", "ascii")
            point_data.appendChild(forces)

            string = self.coords_to_string(x_force, y_force, z_force)
            forceData = doc.createTextNode(string)
            forces.appendChild(forceData)

        # Particle radii
        if len(radii) > 0:
            radiiNode = doc.createElementNS("VTK", "DataArray")
            radiiNode.set_attribute("Name", "radii")
            radiiNode.set_attribute("type", "Float32")
            radiiNode.set_attribute("format", "ascii")
            point_data.appendChild(radiiNode)

            string = self.array_to_string(radii)
            radiiData = doc.createTextNode(string)
            radiiNode.appendChild(radiiData)

        if len(colors) > 0:
            # Particle colors
            colorNode= doc.createElementNS("VTK", "DataArray")
            colorNode.set_attribute("Name", "colors")
            colorNode.set_attribute("type", "Float32")
            colorNode.set_attribute("format", "ascii")
            point_data.appendChild(colorNode)

            string = self.array_to_string(colors)
            color_Data = doc.createTextNode(string)
            colorNode.appendChild(color_Data)

        #### Cell data (dummy) ####
        cell_data = doc.createElementNS("VTK", "CellData")
        piece.appendChild(cell_data)

        # Write to file and exit
        outFile = open(fileName, 'w')
#        xml.dom.ext.PrettyPrint(doc, file)
        doc.writexml(outFile, newl='\n')
        outFile.close()
        self.fileNames.append(fileName)

    def writePVD(self, fileName):
        outFile = open(fileName, 'w')
        import xml.dom.minidom

        pvd = xml.dom.minidom.Document()
        pvd_root = pvd.createElementNS("VTK", "VTKFile")
        pvd_root.set_attribute("type", "Collection")
        pvd_root.set_attribute("version", "0.1")
        pvd_root.set_attribute("byte_order", "LittleEndian")
        pvd.appendChild(pvd_root)

        collection = pvd.createElementNS("VTK", "Collection")
        pvd_root.appendChild(collection)

        for i in range(len(self.fileNames)):
            dataSet = pvd.createElementNS("VTK", "DataSet")
            dataSet.set_attribute("timestep", str(i))
            dataSet.set_attribute("group", "")
            dataSet.set_attribute("part", "0")
            dataSet.set_attribute("file", str(self.fileNames[i]))
            collection.appendChild(dataSet)

        outFile = open(fileName, 'w')
        pvd.writexml(outFile, newl='\n')
        outFile.close()
