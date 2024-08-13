import bpy

bpy.context.scene.render.engine = 'CYCLES'

bpy.ops.object.delete(use_global=False)
#bpy.ops.wm.open_mainfile(filepath="supercell_3d_qc_only_trying_for_python.blend")
bpy.ops.wm.open_mainfile(filepath="supercell_3d_qc_only_straight_cam_for_python.blend")

#bpy.context.scene.cycles.samples = 1024  # Adjust as needed
#bpy.context.scene.render.resolution_x = 8000  # Width
#bpy.context.scene.render.resolution_y = 4000  # Height
bpy.context.scene.render.resolution_percentage = 50  # Full resolution
bpy.context.scene.render.image_settings.file_format='PNG'
bpy.context.scene.render.image_settings.color_depth = '8'  # 8 or 16 bits per channel
bpy.context.scene.render.image_settings.compression = 0  # Set to 15% compression
bpy.context.scene.render.image_settings.color_mode = 'RGBA'

bpy.context.scene.render.use_border = True
bpy.context.scene.render.border_min_x = 0.0
bpy.context.scene.render.border_max_x = 1.0
bpy.context.scene.render.border_min_y = 0.15
bpy.context.scene.render.border_max_y = 0.85

for i in range(5,167,1):

    bpy.ops.object.text_add(location=(0, 0, 0))
    text_obj = bpy.context.object
    text_obj.data.body = "Time (hr): %f"%(i/60.0)
    text_obj.location = (-30, 0, 15)
    text_obj.data.size = 1.0
    text_obj.scale = (2.0, 2.0, 2.0)  # Scale uniformly in all directions
    text_obj.rotation_euler = (90, 0, 0)  # Replace with desired rotation
    text_obj.data.materials.clear()

    if(i<=9):
        qc_filename = 'supercell_3d_00%d0.stl'%i
    elif(i<=99):
        qc_filename = 'supercell_3d_0%d0.stl'%i
    elif(i<=999):
        qc_filename = 'supercell_3d_%d0.stl'%i

    file_loc = './STLFiles/qc/%s'%qc_filename

    if(i<=9):
        img_filename = 'supercell_3d_00%d.jpg'%i
    elif(i<=99):
        img_filename = 'supercell_3d_0%d.jpg'%i
    elif(i<=999):
        img_filename = 'supercell_3d_%d.jpg'%i

    print("filename: ",file_loc);
    #drop=bpy.ops.import_mesh.stl(filepath=file_loc)
    drop=bpy.ops.wm.stl_import(filepath=file_loc)
    mat_drop = bpy.data.materials.get("Drop")
    if mat_drop is None:
        # create material
        mat_drop = bpy.data.materials.new(name="Drop")

    mat_drop.use_nodes=True
    nodes=mat_drop.node_tree.nodes
    for node in nodes:
        nodes.remove(node)

    diffuseshader = nodes.new(type='ShaderNodeBsdfDiffuse')
    diffuseshader.inputs[1].default_value=0.0
    diffuseshader.inputs['Color'].default_value = (0.846, 0.846, 0.846, 1.0)  # RGBA


    #glossyshader = nodes.new(type='ShaderNodeBsdfGlossy')
    #glossyshader.inputs[1].default_value=0.0
    mixshader = nodes.new(type='ShaderNodeMixShader')
    mixshader.inputs['Fac'].default_value = 0.1
    dropoutput = nodes.new(type='ShaderNodeOutputMaterial')

    links = mat_drop.node_tree.links
    link = links.new(diffuseshader.outputs[0], mixshader.inputs[1])
    #link = links.new(glossyshader.outputs[0], mixshader.inputs[2])
    link = links.new(mixshader.outputs[0], dropoutput.inputs[0])

     # Get material
    ob = bpy.context.active_object
    if ob.data.materials:
    # assign to 1st material slot
        ob.data.materials[1] = mat_drop
    else:
    # no slots
        ob.data.materials.append(mat_drop)
    bpy.data.objects[7].scale = (0.001, 0.001, 0.001)

    for area in bpy.context.screen.areas:
        if area.type == 'VIEW_3D':
            area.spaces[0].shading.type = 'RENDERED'

    bpy.context.scene.render.filepath = "./Images/%s"%img_filename
    bpy.ops.render.render(write_still=True)

    bpy.ops.object.select_all(action='DESELECT')
    bpy.data.objects[7].select_set(True)
    bpy.ops.object.delete()
    text_obj.select_set(True)  # Select the object
    bpy.ops.object.delete()  # Delete the selected object
