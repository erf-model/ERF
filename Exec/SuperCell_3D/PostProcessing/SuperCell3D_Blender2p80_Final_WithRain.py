import bpy

bpy.context.scene.render.engine = 'CYCLES'

bpy.ops.object.delete(use_global=False)
#bpy.ops.wm.open_mainfile(filepath="supercell_3d_qc_only_trying_for_python.blend")
#bpy.ops.wm.open_mainfile(filepath="supercell_3d_qc_only_straight_cam_for_python.blend")
bpy.ops.wm.open_mainfile(filepath="supercell_3d_qc_only_8k_for_python_black.blend")

#bpy.context.scene.cycles.samples = 1024  # Adjust as needed
bpy.context.scene.render.resolution_x = 7680  # Width
bpy.context.scene.render.resolution_y = 4320  # Height
bpy.context.scene.render.resolution_percentage = 50  # Full resolution   
bpy.context.scene.render.image_settings.file_format='PNG'
bpy.context.scene.render.image_settings.color_depth = '8'  # 8 or 16 bits per channel
bpy.context.scene.render.image_settings.compression = 0  # Set to 15% compression
bpy.context.scene.render.image_settings.color_mode = 'RGBA'

bpy.context.scene.render.use_border = True
#bpy.context.scene.render.border_min_x = 0.0
#bpy.context.scene.render.border_max_x = 0.0
bpy.context.scene.render.border_min_y = 0.0
bpy.context.scene.render.border_max_y = 0.6

camera = bpy.context.scene.camera

# Set the Clip Start and Clip End values
camera.data.clip_end = 150.0  # Example: start clipping at 0.1 Blender units

for i in range(5,167,1):
	
	bpy.ops.object.text_add(location=(0, 0, 0))
	text_obj1 = bpy.context.object
	text_obj1.data.body = "Time (hr): %0.2f"%(i/60.0)
	text_obj1.location = (-32, 0, 17)
	text_obj1.data.size = 1.0
	text_obj1.scale = (2.0, 2.0, 2.0)  # Scale uniformly in all directions
	text_obj1.rotation_euler = (90, 0, 0)  # Replace with desired rotation

	# Create a new material
	text_material = bpy.data.materials.new(name="TextMaterial")

	# Enable 'Use Nodes' to allow color changes
	text_material.use_nodes = True

	# Access the material's nodes and set the color
	bsdf = text_material.node_tree.nodes["Principled BSDF"]
	bsdf.inputs['Base Color'].default_value = (1.0, 1.0, 1.0, 1.0)  # Red color (RGBA)

	# Assign the material to the text object
	if text_obj1.data.materials:
		text_obj1.data.materials[0] = text_material
	else:
		text_obj1.data.materials.append(text_material)

	#text_obj.data.materials.clear()

	bpy.ops.object.text_add(location=(0, 0, 0))
	text_obj = bpy.context.object
	text_obj.data.body = "Scale"
	text_obj.location = (10.0, 0, 17.25)
	text_obj.data.size = 0.75
	text_obj.scale = (2.0, 2.0, 2.0)  # Scale uniformly in all directions
	text_obj.rotation_euler = (90, 0, 0)  # Replace with desired rotation
	#text_obj.data.materials.clear()

	# Create a new material
	text_material = bpy.data.materials.new(name="TextMaterial")

	# Enable 'Use Nodes' to allow color changes
	text_material.use_nodes = True

	# Access the material's nodes and set the color
	bsdf = text_material.node_tree.nodes["Principled BSDF"]
	bsdf.inputs['Base Color'].default_value = (1.0, 1.0, 1.0, 1.0)  # Red color (RGBA)

	# Assign the material to the text object
	if text_obj.data.materials:
		text_obj.data.materials[0] = text_material
	else:
		text_obj.data.materials.append(text_material)


	bpy.ops.object.text_add(location=(0, 0, 0))
	text_obj = bpy.context.object
	text_obj.data.body = "__"
	text_obj.location = (13.5, 0, 18.4)
	text_obj.data.size = 5.0
	text_obj.scale = (2.0, 2.0, 2.0)  # Scale uniformly in all directions
	text_obj.rotation_euler = (90, 0, 0)  # Replace with desired rotation
	#text_obj.data.materials.clear()

	# Create a new material
	text_material = bpy.data.materials.new(name="TextMaterial")

	# Enable 'Use Nodes' to allow color changes
	text_material.use_nodes = True

	# Access the material's nodes and set the color
	bsdf = text_material.node_tree.nodes["Principled BSDF"]
	bsdf.inputs['Base Color'].default_value = (1.0, 1.0, 1.0, 1.0)  # Red color (RGBA)

	# Assign the material to the text object
	if text_obj.data.materials:
		text_obj.data.materials[0] = text_material
	else:
		text_obj.data.materials.append(text_material)


	bpy.ops.object.text_add(location=(0, 0, 0))
	text_obj = bpy.context.object
	text_obj.data.body = "10 km"
	text_obj.location = (16.7, 0, 18.0)
	text_obj.data.size = 0.75
	text_obj.scale = (2.0, 2.0, 2.0)  # Scale uniformly in all directions
	text_obj.rotation_euler = (90, 0, 0)  # Replace with desired rotation
	text_obj.data.materials.clear()

		# Create a new material
	text_material = bpy.data.materials.new(name="TextMaterial")

	# Enable 'Use Nodes' to allow color changes
	text_material.use_nodes = True

	# Access the material's nodes and set the color
	bsdf = text_material.node_tree.nodes["Principled BSDF"]
	bsdf.inputs['Base Color'].default_value = (1.0, 1.0, 1.0, 1.0)  # Red color (RGBA)

	# Assign the material to the text object
	if text_obj.data.materials:
		text_obj.data.materials[0] = text_material
	else:
		text_obj.data.materials.append(text_material)



	if(i<=9):
		qc_filename = 'supercell_3d_00%d0.stl'%i
	elif(i<=99):
		qc_filename = 'supercell_3d_0%d0.stl'%i
	elif(i<=999):
		qc_filename = 'supercell_3d_%d0.stl'%i

	file_loc = './STLFiles/qc/%s'%qc_filename 
	file_loc_qr = './STLFiles/qr/%s'%qc_filename 

	if(i<=9):
		img_filename = 'supercell_3d_00%d.jpg'%i
	elif(i<=99):
		img_filename = 'supercell_3d_0%d.jpg'%i
	elif(i<=999):
		img_filename = 'supercell_3d_%d.jpg'%i
	
	print("filename: ",file_loc);
	drop=bpy.ops.wm.stl_import(filepath=file_loc)
	obj1 = bpy.context.selected_objects[0]  # Get the imported object
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
	diffuseshader.inputs['Color'].default_value = (0.0, 0.0, 0.906, 1.0)  # RGBA

	#glossyshader = nodes.new(type='ShaderNodeBsdfGlossy')
	#glossyshader.inputs[1].default_value=0.0
	mixshader = nodes.new(type='ShaderNodeMixShader')
	mixshader.inputs['Fac'].default_value = 0.1
	dropoutput = nodes.new(type='ShaderNodeOutputMaterial')

	principledshader = nodes.new(type='ShaderNodeBsdfPrincipled')
	principledshader.inputs['Alpha'].default_value = 1.0
	#principledshader.inputs['Base Color'].default_value = 1.0  # Red color (RGBA)

	links = mat_drop.node_tree.links
	link = links.new(diffuseshader.outputs[0], mixshader.inputs[1])
	#link = links.new(glossyshader.outputs[0], mixshader.inputs[2])
	link = links.new(mixshader.outputs[0], principledshader.inputs[0])
	link = links.new(principledshader.outputs[0], dropoutput.inputs[0])

	# Rain water
	if(i >=10):
		rain=bpy.ops.wm.stl_import(filepath=file_loc_qr)
		obj2 = bpy.context.selected_objects[0]  # Get the imported object
		mat_rain = bpy.data.materials.get("Rain")
	
		if mat_rain is None:
			# create material
			mat_rain = bpy.data.materials.new(name="Rain")

		mat_rain.use_nodes=True
		nodes=mat_rain.node_tree.nodes
		for node in nodes:
			nodes.remove(node)

		glossyshader = nodes.new(type='ShaderNodeBsdfGlossy')
		glossyshader.inputs[1].default_value=0.0
		glossyshader.inputs['Color'].default_value = (0.17, 0.481, 0.8, 1.0)  # RGBA

		refractionshader = nodes.new(type='ShaderNodeBsdfRefraction')
		refractionshader.inputs[1].default_value=0.0
		refractionshader.inputs['Color'].default_value = (0.089, 0.299, 0.8, 1.0)  # RGBA

		mixshader = nodes.new(type='ShaderNodeMixShader')
		mixshader.inputs['Fac'].default_value = 0.5
		dropoutput = nodes.new(type='ShaderNodeOutputMaterial')

		links = mat_rain.node_tree.links
		link = links.new(refractionshader.outputs[0], mixshader.inputs[1])
		link = links.new(glossyshader.outputs[0], mixshader.inputs[2])
		link = links.new(mixshader.outputs[0], dropoutput.inputs[0])

	
	 # Get material
	#ob = bpy.context.active_object
	if obj1.data.materials:
	# assign to 1st material slot
		obj1.data.materials[0] = mat_drop
	else:
	# no slots
		obj1.data.materials.append(mat_drop)

	if(i>=10):

		if obj2.data.materials:
		# assign to 1st material slot
			obj2.data.materials[0] = mat_rain
		else:
		# no slots
			obj2.data.materials.append(mat_rain)

	obj1.scale = (0.001, 0.001, 0.001)
	if(i>=10):
		obj2.scale = (0.001, 0.001, 0.001)

	for area in bpy.context.screen.areas:
		if area.type == 'VIEW_3D':
			area.spaces[0].shading.type = 'RENDERED'
	
	bpy.context.scene.render.filepath = "./Images/%s"%img_filename
	bpy.ops.render.render(write_still=True)

	# Remove the materials after rendering
	bpy.data.materials.remove(mat_drop)
	if(i>=10):
		bpy.data.materials.remove(mat_rain)

	# Optionally, delete the imported objects
	bpy.data.objects.remove(obj1)
	if(i>=10):
		bpy.data.objects.remove(obj2)
	
	#bpy.ops.object.select_all(action='DESELECT')
	#bpy.data.objects[7].select_set(True)
	#bpy.data.objects[8].select_set(True)
	#bpy.ops.object.delete()

	# After you're done with the text object, delete it
	bpy.context.view_layer.objects.active = text_obj1  # Set as active
	text_obj1.select_set(True)  # Select the object
	bpy.ops.object.delete()  # Delete the selected object

