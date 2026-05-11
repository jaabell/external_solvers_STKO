import PyMpc.Units as u
from PyMpc import *
from mpc_utils_html import *
import opensees.utils.tcl_input as tclin

def _err(msg):
	return 'Error in "stressControl" :\n{}'.format(msg)

def _geta(xobj, name):
	a = xobj.getAttribute(name)
	if a is None:
		raise Exception(_err('cannot find "{}" attribute'.format(name)))
	return a

def makeXObjectMetaData():

	def mka(type, name, group, descr):
		a = MpcAttributeMetaData()
		a.type = type
		a.name = name
		a.group = group
		a.description = (
			html_par(html_begin()) +
			html_par(html_boldtext(name)+'<br/>') +
			html_par(descr) +
			html_end()
			)
		return a

	sset = mka(MpcAttributeType.IndexVector, 'SelectionSets', 'Selection',
		'All selection sets whose elements are to be used for stress control.')
	sset.indexSource.type = MpcAttributeIndexSourceType.SelectionSet

	auto_gen = mka(MpcAttributeType.Boolean, 'Include Auto-Generated Elements', 'Selection',
		('Some automations in STKO will automatically generate extra elements not visible in STKO (see for example the HingedBeam).<br/>'
		'When this flag is True (Default), the auto-generated elements will be considered.'))
	auto_gen.setDefault(True)

	duration = mka(MpcAttributeType.Real, 'Total Duration', 'Time',
		'Total analysis duration. The target stress increments are reached when the pseudo-time reaches this value.')
	duration.setDefault(1.0)

	nsteps = mka(MpcAttributeType.Integer, 'Number of Steps', 'Time',
		'Number of analysis steps over which the stress ramp is applied.')
	nsteps.setDefault(10)

	target_xx = mka(MpcAttributeType.Real, 'Target Stress Increment XX', 'Targets',
		'Final increment of commitStressXX reached when pseudo-time = Total Duration.')
	target_xx.setDefault(0.0)
	target_yy = mka(MpcAttributeType.Real, 'Target Stress Increment YY', 'Targets',
		'Final increment of commitStressYY reached when pseudo-time = Total Duration.')
	target_yy.setDefault(0.0)
	target_zz = mka(MpcAttributeType.Real, 'Target Stress Increment ZZ', 'Targets',
		'Final increment of commitStressZZ reached when pseudo-time = Total Duration.')
	target_zz.setDefault(0.0)
	target_xy = mka(MpcAttributeType.Real, 'Target Stress Increment XY', 'Targets',
		'Final increment of commitStressXY reached when pseudo-time = Total Duration.')
	target_xy.setDefault(0.0)
	target_yz = mka(MpcAttributeType.Real, 'Target Stress Increment YZ', 'Targets',
		'Final increment of commitStressYZ reached when pseudo-time = Total Duration.')
	target_yz.setDefault(0.0)
	target_xz = mka(MpcAttributeType.Real, 'Target Stress Increment XZ', 'Targets',
		'Final increment of commitStressXZ reached when pseudo-time = Total Duration.')
	target_xz.setDefault(0.0)

	xom = MpcXObjectMetaData()
	xom.name = 'stressControl'
	xom.addAttribute(sset)
	xom.addAttribute(auto_gen)
	xom.addAttribute(duration)
	xom.addAttribute(nsteps)
	xom.addAttribute(target_xx)
	xom.addAttribute(target_yy)
	xom.addAttribute(target_zz)
	xom.addAttribute(target_xy)
	xom.addAttribute(target_yz)
	xom.addAttribute(target_xz)
	return xom

def writeTcl(pinfo):

	from io import StringIO

	doc = App.caeDocument()
	xobj = pinfo.analysis_step.XObject

	sset = _geta(xobj, 'SelectionSets').indexVector
	auto_gen = _geta(xobj, 'Include Auto-Generated Elements').boolean
	total_duration = _geta(xobj, 'Total Duration').real
	num_steps = _geta(xobj, 'Number of Steps').integer

	targets = [
		('XX', _geta(xobj, 'Target Stress Increment XX').real),
		('YY', _geta(xobj, 'Target Stress Increment YY').real),
		('ZZ', _geta(xobj, 'Target Stress Increment ZZ').real),
		('XY', _geta(xobj, 'Target Stress Increment XY').real),
		('XZ', _geta(xobj, 'Target Stress Increment XZ').real),
		('YZ', _geta(xobj, 'Target Stress Increment YZ').real),
	]

	# filter out zero targets
	targets = [(comp, val) for comp, val in targets if val != 0.0]

	if len(targets) == 0:
		return

	# collect elements from selection sets
	parameter_map_elem = []
	auto_gen_elements_pid_map = {}

	def add_elements_to_param(source_elements):
		for elem in source_elements:
			parameter_map_elem.append(elem.id)
			if auto_gen:
				aux_elements = pinfo.auto_generated_element_data_map.get(elem.id, None)
				source_pid = 0
				if pinfo.process_count > 1:
					source_pid = doc.mesh.partitionData.elementPartition(elem.id)
				if aux_elements is not None:
					for aux_ele_id in aux_elements.elements:
						parameter_map_elem.append(aux_ele_id)
						if pinfo.process_count > 1:
							auto_gen_elements_pid_map[aux_ele_id] = source_pid

	for selection_set_id in sset:
		if not selection_set_id in doc.selectionSets: continue
		selection_set = doc.selectionSets[selection_set_id]
		for geometry_id, geometry_subset in selection_set.geometries.items():
			mesh_of_geom = doc.mesh.meshedGeometries[geometry_id]
			for domain_id in geometry_subset.edges:
				domain = mesh_of_geom.edges[domain_id]
				add_elements_to_param(domain.elements)
			for domain_id in geometry_subset.faces:
				domain = mesh_of_geom.faces[domain_id]
				add_elements_to_param(domain.elements)
			for domain_id in geometry_subset.solids:
				domain = mesh_of_geom.solids[domain_id]
				add_elements_to_param(domain.elements)
		for interaction_id in selection_set.interactions:
			mesh_of_inter = doc.mesh.meshedInteractions[interaction_id]
			add_elements_to_param(mesh_of_inter.elements)

	# partition elements by process id
	pid_element_map = {}
	if pinfo.process_count > 1:
		for ele_id in parameter_map_elem:
			if ele_id in auto_gen_elements_pid_map:
				pid = auto_gen_elements_pid_map[ele_id]
			else:
				pid = doc.mesh.partitionData.elementPartition(ele_id)
			pid_values = pid_element_map.get(pid, None)
			if pid_values is None:
				pid_values = []
				pid_element_map[pid] = pid_values
			pid_values.append(ele_id)
	else:
		pid_element_map[0] = parameter_map_elem

	# generate unique parameter tags for each component
	comp_param_tags = []
	for comp, val in targets:
		tag = pinfo.tag_parameter
		pinfo.tag_parameter += 1
		comp_param_tags.append(tag)

	step_id = pinfo.analysis_step.id
	proc_name = '_stressCtrl_{}'.format(step_id)
	var_prefix = '_stressCtrl_{}'.format(step_id)

	# write header comment
	pinfo.out_file.write('\n{}# stressControl (step {})\n'.format(pinfo.indent, step_id))

	# build element list in TCL (one list per process)
	if pinfo.process_count > 1:
		for pid, elements in pid_element_map.items():
			if len(elements) > 0:
				pinfo.out_file.write('{}if {{$STKO_VAR_process_id == {}}} {{\n'.format(pinfo.indent, pid))
				ele_list = ' '.join(str(e) for e in elements)
				pinfo.out_file.write('{}{}set {}_elems_{} {{ {} }}\n'.format(
					pinfo.indent, pinfo.tabIndent, var_prefix, pid, ele_list))
				pinfo.out_file.write('{}}}\n'.format(pinfo.indent))
	else:
		elements = pid_element_map[0]
		ele_list = ' '.join(str(e) for e in elements)
		pinfo.out_file.write('{}set {}_elems {{ {} }}\n'.format(pinfo.indent, var_prefix, ele_list))

	# write parameter and addToParameter commands using foreach
	for (comp, val), tag in zip(targets, comp_param_tags):
		param_name = 'commitStressIncrement' + comp
		pinfo.out_file.write('{}parameter {}\n'.format(pinfo.indent, tag))
		if pinfo.process_count > 1:
			for pid, _ in pid_element_map.items():
				pinfo.out_file.write('{}if {{$STKO_VAR_process_id == {}}} {{\n'.format(pinfo.indent, pid))
				pinfo.out_file.write('{}{}foreach _stressCtrl_elem ${}_elems_{} {{\n'.format(
					pinfo.indent, pinfo.tabIndent, var_prefix, pid))
				pinfo.out_file.write('{}{}{}addToParameter {} element $_stressCtrl_elem "{}"\n'.format(
					pinfo.indent, pinfo.tabIndent, pinfo.tabIndent, tag, param_name))
				pinfo.out_file.write('{}{}}}\n'.format(pinfo.indent, pinfo.tabIndent))
				pinfo.out_file.write('{}}}\n'.format(pinfo.indent))
		else:
			pinfo.out_file.write('{}foreach _stressCtrl_elem ${}_elems {{\n'.format(pinfo.indent, var_prefix))
			pinfo.out_file.write('{}{}addToParameter {} element $_stressCtrl_elem "{}"\n'.format(
				pinfo.indent, pinfo.tabIndent, tag, param_name))
			pinfo.out_file.write('{}}}\n'.format(pinfo.indent))

	# --- write the OnBeforeAnalyze proc ---
	pinfo.out_file.write('\n{}# stressControl proc (step {})\n'.format(pinfo.indent, step_id))
	pinfo.out_file.write('{}proc {} {{}} {{\n'.format(pinfo.indent, proc_name))
	ind = pinfo.indent + pinfo.tabIndent

	# tracking array variable name
	var_prefix = '_stressCtrl_{}'.format(step_id)

	# init tracking on first call
	pinfo.out_file.write('{}global {}\n'.format(ind, var_prefix))
	pinfo.out_file.write('{}if {{![info exists {}(count)]}} {{\n'.format(ind, var_prefix))
	pinfo.out_file.write('{}{}set {}(count) 0\n'.format(ind, pinfo.tabIndent, var_prefix))
	for comp, val in targets:
		pinfo.out_file.write('{}{}set {}({}) 0.0\n'.format(ind, pinfo.tabIndent, var_prefix, comp))
	pinfo.out_file.write('{}}}\n'.format(ind))

	# increment step counter and compute factor
	pinfo.out_file.write('{}set {}(count) [expr {{${}(count)}} + 1]\n'.format(ind, var_prefix, var_prefix))
	pinfo.out_file.write('{}set _stressCtrl_factor [expr {{${}(count) / {}.0}}]\n'.format(ind, var_prefix, num_steps))
	pinfo.out_file.write('{}if {{$_stressCtrl_factor > 1.0}} {{ set _stressCtrl_factor 1.0 }}\n'.format(ind))

	# for each component
	for (comp, val), tag in zip(targets, comp_param_tags):
		pinfo.out_file.write('{}set _stressCtrl_current [expr {} * $_stressCtrl_factor]\n'.format(ind, val))
		pinfo.out_file.write('{}set _stressCtrl_incr [expr $_stressCtrl_current - ${}({})]\n'.format(ind, var_prefix, comp))
		pinfo.out_file.write('{}updateParameter {} $_stressCtrl_incr\n'.format(ind, tag))
		pinfo.out_file.write('{}set {}({}) $_stressCtrl_current\n'.format(ind, var_prefix, comp))

	pinfo.out_file.write('{}}}\n'.format(pinfo.indent))  # end proc

	# register in OnBeforeAnalyze
	pinfo.out_file.write('{}lappend STKO_VAR_OnBeforeAnalyze_CustomFunctions {}\n'.format(
		pinfo.indent, proc_name))
