import PyMpc.Units as u
from PyMpc import *
from mpc_utils_html import *
import opensees.utils.tcl_input as tclin

def makeXObjectMetaData():
	# -- Copy and reuse all attributes as-is --
	at_thick = MpcAttributeMetaData()
	at_thick.type = MpcAttributeType.QuantityScalar
	at_thick.name = 'thick'
	at_thick.group = 'Group'
	at_thick.description = (
		html_par(html_begin()) +
		html_par(html_boldtext('thick') + '<br/>') +
		html_par('element thickness') +
		html_par(html_href('https://opensees.berkeley.edu/wiki/index.php/tri6n_Element', 'tri6n Element') + '<br/>') +
		html_end()
	)
	at_thick.dimension = u.L

	at_type = MpcAttributeMetaData()
	at_type.type = MpcAttributeType.String
	at_type.name = 'type'
	at_type.group = 'Group'
	at_type.description = (
		html_par(html_begin()) +
		html_par(html_boldtext('type') + '<br/>') +
		html_par('string representing material behavior. The type parameter can be either "PlaneStrain" or "PlaneStress."') +
		html_par(html_href('https://opensees.berkeley.edu/wiki/index.php/tri6n_Element', 'tri6n Element') + '<br/>') +
		html_end()
	)
	at_type.sourceType = MpcAttributeSourceType.List
	at_type.setSourceList(['PlaneStrain', 'PlaneStress'])
	at_type.setDefault('PlaneStrain')

	at_Optional = MpcAttributeMetaData()
	at_Optional.type = MpcAttributeType.Boolean
	at_Optional.name = 'Optional'
	at_Optional.group = 'Group'
	at_Optional.description = (
		html_par(html_begin()) +
		html_par(html_boldtext('Optional') + '<br/>') +
		html_par('to activate pressure, rho, b1 and b2') +
		html_end()
	)

	at_pressure = MpcAttributeMetaData()
	at_pressure.type = MpcAttributeType.QuantityScalar
	at_pressure.name = 'pressure'
	at_pressure.group = 'Optional parameters'
	at_pressure.description = (
		html_par(html_begin()) +
		html_par(html_boldtext('pressure') + '<br/>') +
		html_par('surface pressure (optional, default = 0.0)') +
		html_end()
	)
	at_pressure.setDefault(0.0)
	at_pressure.dimension = u.F/u.L**2

	at_rho = MpcAttributeMetaData()
	at_rho.type = MpcAttributeType.QuantityScalar
	at_rho.name = 'rho'
	at_rho.group = 'Optional parameters'
	at_rho.description = (
		html_par(html_begin()) +
		html_par(html_boldtext('rho') + '<br/>') +
		html_par('element mass density (optional, default = 0.0)') +
		html_end()
	)
	at_rho.setDefault(0.0)

	at_b1 = MpcAttributeMetaData()
	at_b1.type = MpcAttributeType.QuantityScalar
	at_b1.name = 'b1'
	at_b1.group = 'Optional parameters'
	at_b1.description = (
		html_par(html_begin()) +
		html_par(html_boldtext('b1') + '<br/>') +
		html_par('body force in x-direction (optional, default = 0.0)') +
		html_end()
	)
	at_b1.setDefault(0.0)
	at_b1.dimension = u.F

	at_b2 = MpcAttributeMetaData()
	at_b2.type = MpcAttributeType.QuantityScalar
	at_b2.name = 'b2'
	at_b2.group = 'Optional parameters'
	at_b2.description = (
		html_par(html_begin()) +
		html_par(html_boldtext('b2') + '<br/>') +
		html_par('body force in y-direction (optional, default = 0.0)') +
		html_end()
	)
	at_b2.setDefault(0.0)
	at_b2.dimension = u.F

	xom = MpcXObjectMetaData()
	xom.name = 'tri6n'
	xom.addAttribute(at_thick)
	xom.addAttribute(at_type)
	xom.addAttribute(at_Optional)
	xom.addAttribute(at_pressure)
	xom.addAttribute(at_rho)
	xom.addAttribute(at_b1)
	xom.addAttribute(at_b2)

	xom.setVisibilityDependency(at_Optional, at_pressure)
	xom.setVisibilityDependency(at_Optional, at_rho)
	xom.setVisibilityDependency(at_Optional, at_b1)
	xom.setVisibilityDependency(at_Optional, at_b2)

	return xom

def getNodalSpatialDim(xobj, xobj_phys_prop):
	return [(2,2)] * 6

def writeTcl(pinfo):
	# element tri6n $tag $n1 $n2 $n3 $n4 $n5 $n6 $thick $type $matTag <$pressure $rho $b1 $b2>
	elem = pinfo.elem
	phys_prop = pinfo.phys_prop
	elem_prop = pinfo.elem_prop

	tag = elem.id
	matTag = phys_prop.id
	xobj = elem_prop.XObject

	ClassName = xobj.name
	if pinfo.currentDescription != ClassName:
		pinfo.out_file.write('\n{}# {} {}\n'.format(pinfo.indent, xobj.Xnamespace, ClassName))
		pinfo.currentDescription = ClassName

	namePh = phys_prop.XObject.Xnamespace
	if not namePh.startswith('materials.nD'):
		raise Exception('Error: physical property must be "materials.nD" and not: "{}"'.format(namePh))

	# mandatory
	thick = xobj.getAttribute('thick').quantityScalar
	type = xobj.getAttribute('type').string

	# validate geometry
	if (elem.geometryFamilyType() != MpcElementGeometryFamilyType.Triangle) or len(elem.nodes) != 6:
		raise Exception('Error: element must be a 6-node triangle')

	node_ids = [node.id for node in elem.nodes]
	node_str = ' '.join(str(nid) for nid in node_ids)

	# optional parameters
	sopt = ''
	if xobj.getAttribute('Optional').boolean:
		pressure = xobj.getAttribute('pressure').quantityScalar.value
		rho = xobj.getAttribute('rho').quantityScalar.value
		b1 = xobj.getAttribute('b1').quantityScalar.value
		b2 = xobj.getAttribute('b2').quantityScalar.value
		sopt = ' {} {} {} {}'.format(pressure, rho, b1, b2)

	# write tcl
	pinfo.out_file.write('{}element tri6n {} {} {} {} {}\n'.format(
		pinfo.indent, tag, node_str, thick.value, type, matTag) + sopt + '\n')
