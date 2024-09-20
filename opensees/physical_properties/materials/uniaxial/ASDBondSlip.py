# enable default 3D tester for this module
from opensees.physical_properties.utils.tester.EnableTester1D import *
from opensees.utils.override_utils import get_function_from_module

import PyMpc.Units as u
from PyMpc import *
from mpc_utils_html import *
import opensees.utils.tcl_input as tclin
import math

def _err(msg):
	return 'Error in ASDBondSlip: {}'.format(msg)

def _geta(xobj, name):
	at = xobj.getAttribute(name)
	if at is None:
		raise Exception(_err('cannot find "{}" attribute'.format(name)))
	return at

class _globals:
	L_units = {
		'mm' : 1.0,
		'cm' : 10.0,
		'dm' : 100.0,
		'm'  : 1000.0,
		'in' : 25.4,
		'ft' : 304.8,
		}
	F_units = {
		'N': 1.0,
		'daN': 10.0,
		'kN' : 1000.0,
		'lbf' : 4.4482216152605,
		'kip' : 4448.2216152605,
		}
	concrete_format = ('{0}uniaxialMaterial ASDConcrete1D {1} {2} \\\n'
		'{0}\t-Te {4} \\\n'
		'{0}\t-Ts {5} \\\n'
		'{0}\t-Td {6} \\\n'
		'{0}\t-Ce {7} \\\n'
		'{0}\t-Cs {8} \\\n'
		'{0}\t-Cd {9} \\\n'
		'{0}\t -eta {3}{10}{11}\n')

class _mc2020:
	alpha = 0.4
	tolerance = 1.0e-4
	def _discretize(ymax,s1):
		# dicretize the first part of the equation in 4 points
		Y = [0.0, 0.3*ymax, 0.6*ymax, 0.8*ymax]
		X = [(ti/ymax)**(1/_mc2020.alpha)*s1 for ti in Y]
		return (X,Y)
	def _make_po_good(fc, cclear):
		tmax = 2.5*math.sqrt(fc)
		s1 = 1.0
		s2 = 2.0
		s3 = cclear if cclear > s2 else s2+0.01
		tu = 0.4*tmax
		X,Y = _mc2020._discretize(tmax,s1)
		X.extend([s1, s2, s3])
		Y.extend([tmax, tmax, tu])
		return (X,Y)
	def make(pull_out, good_bond, confined, fc, cclear):
		# obtain discrete points of the total backbone curve
		X,Y = _mc2020._make_po_good(fc, cclear)
		# according to _discretize function, the peak force is at the 5-th entry
		tpeak = Y[4]
		xpeak = X[4]
		tu = Y[-1]
		stress_null = tpeak*_mc2020.tolerance
		# now create the pinching part, and the frictional residual parts
		# so that their sum is equal to the total backbone
		Y1 = [yi-tu*(xi/xpeak) if xi < xpeak else yi-tu for xi,yi in zip(X,Y)]
		Y2 = [tu*(xi/xpeak) if xi < xpeak else tu for xi,yi in zip(X,Y)]
		# according to the mc2020, the unloading modulus is secant up to 0.6*tmax (3rd point in _discretize)
		xdam = X[2]
		D = [1.0 if xi <= xdam else 0.0 for xi in X]
		# done
		return (X, Y1, Y2, D, stress_null)

def onAttributeChanged(editor, xobj, attribute_name):
	if attribute_name == 'Failure Mode':
		is_splitting = _geta(xobj, attribute_name).string == 'Splitting(SP)'
		_geta(xobj, 'Confinement').visible = is_splitting
	return None

def makeXObjectMetaData():
	
	def mka(name, group, descr, atype, adim = None, dval = None):
		a = MpcAttributeMetaData()
		a.type = atype
		a.name = name
		a.group = group
		a.description = (
			html_par(html_begin()) +
			html_par(html_boldtext(name)+'<br/>') + 
			html_par(descr) +
			#html_par(html_href('https://opensees.github.io/OpenSeesDocumentation/user/manual/material/uniaxialMaterials/ASDConcrete1D.html','ASDConcrete1D')+'<br/>') +
			html_end()
			)
		if adim is not None:
			a.dimension = adim
		if dval is not None:
			a.setDefault(dval)
		return a
	
	fail_type = mka('Failure Mode', 'Bond-Slip', 
		'''<p>"<strong>Pull-out failure</strong>" is valid for well-confined concrete (concrete cover&nbsp;&ge; 5&Oslash;, clear spacing between bars &ge; 10&Oslash;), or suitable confining reinforcement.</p>
<p><span style="color: #3366ff;"><em>Ref: fib Model Code 2020; 20.5.1.1 (pag. 308, Table 20.5-1)</em></span></p>''', 
		MpcAttributeType.String, dval='Splitting(SP)')
	fail_type.sourceType = MpcAttributeSourceType.List
	fail_type.setSourceList(['Splitting(SP)', 'Pull-out(PO)'])
	
	bond_type = mka('Bond conditions', 'Bond-Slip',
'''<p>"<strong>Good</strong>" bond conditions are defined as:</p>
<ul>
<li>bars with an inclination less than 45&deg; to the horizontal which are up to 300 mm from the bottom of the formwork&nbsp;or at least 300 mm below the free surface during concreting.</li>
</ul>
<p>Otherwise "<strong>Poor</strong>" bond conditions should be assumed.</p>
<p><span style="color: #3366ff;"><em>Ref: fib Model Code 2020; 20.2.2.1 (pag. 280)</em></span></p>''',
		MpcAttributeType.String, dval='Poor')
	bond_type.sourceType = MpcAttributeSourceType.List
	bond_type.setSourceList(['Good', 'Poor'])
	
	conf_type = mka('Confinement', 'Bond-Slip', 
		'''<p>"<strong>Unconfined</strong>" or "<strong>Sitrrups</strong>" options available only for "<strong>Splitting (SP)</strong>" failure mode.</p>
<p><span style="color: #3366ff;"><em>Ref: fib Model Code 2020; 20.5.1.1 (pag. 308, Table 20.5-1)</em></span></p>''', 
		MpcAttributeType.String, dval='Unconfined')
	conf_type.sourceType = MpcAttributeSourceType.List
	conf_type.setSourceList(['Unconfined', 'Stirrups'])
	
	fc = mka('fc', 'Bond-Slip', 
		'''<p>Concrete compressive strength</p>
<p><span style="color: #3366ff;"><em>Ref: fib Model Code 2020; 20.5.1.1 (pag. 308, Table 20.5-1)</em></span></p>''', 
		MpcAttributeType.QuantityScalar, adim=u.F/u.L**2, dval=0.0)
	
	cclear = mka('cclear', 'Bond-Slip', 
		'''<p>The clear distance between ribs</p>
<p><span style="color: #3366ff;"><em>Ref: fib Model Code 2020; 20.5.1.1 (pag. 308, Table 20.5-1)</em></span></p>''', 
		MpcAttributeType.QuantityScalar, adim=u.L, dval=0.0)
	
	Lunit = mka('L. unit', 'Units', 'Unit of measurement used for Length', MpcAttributeType.String, dval="mm")
	Lunit.sourceType = MpcAttributeSourceType.List
	Lunit.setSourceList(list(_globals.L_units.keys()))
	
	Funit = mka('F. unit', 'Units', 'Unit of measurement used for Force', MpcAttributeType.String, dval="N")
	Funit.sourceType = MpcAttributeSourceType.List
	Funit.setSourceList(list(_globals.F_units.keys()))
	
	algo = mka('Integration', 'Misc', '''<p>The integration algorithm</p>
<ul>
<li><strong>Implicit</strong>: A standard Backward-Euler integration scheme.</li>
<li><strong>IMPL-EX</strong>:&nbsp;A mixed IMPLicit-EXplicit integration scheme. The resulting response is&nbsp;step-wise linear with a positive-definite tangent stiffness matrix due to the explicit extrapolation of the internal variables.&nbsp;However, the time-step should be smaller than the one used for an implicit scheme.</li>
</ul>
<p><span style="color: #3366ff;"><em>Ref: Oliver, J., Huespe, A. E., &amp; Cante, J. C. (2008). &ldquo;An implicit/explicit integration scheme to increase computability of non-linear material and contact/friction problems&rdquo; Computer Methods in Applied Mechanics and Engineering, 197(21-24), 1865-1889</em></span></p>
<p><span style="color: #3366ff;"><em><a href="https://core.ac.uk/download/pdf/325948712.pdf">Link to article</a></em></span></p>''', 
		MpcAttributeType.String, dval='Implicit')
	algo.sourceType = MpcAttributeSourceType.List
	algo.setSourceList(['Implicit', 'IMPL-EX'])
	
	ctype = mka('Constitutive Tensor Type', 'Misc', '''<p>Constitutive tensor type</p>
<ul>
<li><strong>Tangent</strong>: The algorithmic tangent tensor.</li>
<li><strong>Secant</strong>: The secant tensor.</li>
</ul>
<p>&nbsp;</p>''', 
		MpcAttributeType.String, dval='Secant')
	ctype.sourceType = MpcAttributeSourceType.List
	ctype.setSourceList(['Tangent', 'Secant'])
	
	eta  = mka('eta', 'Misc', 'Viscosity parameter', MpcAttributeType.Real, dval=0.0)
	
	xom = MpcXObjectMetaData()
	xom.name = 'ASDBondSlip'
	xom.Xgroup = 'ASDEASoftware'
	
	xom.addAttribute(fail_type)
	xom.addAttribute(bond_type)
	xom.addAttribute(conf_type)
	xom.addAttribute(fc)
	xom.addAttribute(cclear)
	xom.addAttribute(Lunit)
	xom.addAttribute(Funit)
	xom.addAttribute(eta)
	xom.addAttribute(algo)
	xom.addAttribute(ctype)
	
	return xom

def writeTcl(pinfo):
	
	# xobject and material tag
	xobj = pinfo.phys_prop.XObject
	tag = xobj.parent.componentId
	
	# get basic parameters
	pull_out = _geta(xobj, 'Failure Mode').string == 'Pull-out(PO)'
	good_bond = _geta(xobj, 'Bond conditions').string == 'Good'
	confined = _geta(xobj, 'Confinement').string == 'Stirrups'
	fc = _geta(xobj, 'fc').quantityScalar.value
	cclear = _geta(xobj, 'cclear').quantityScalar.value
	eta = _geta(xobj, 'eta').real
	tangent = _geta(xobj, 'Constitutive Tensor Type').string == 'Tangent'
	implex = _geta(xobj, 'Integration').string == 'IMPL-EX'
	
	# checks
	if fc <= 0.0:
		raise Exception(_err('fc should be strictly positive (>=0)'))
	if cclear <= 0.0:
		raise Exception(_err('cclear should be strictly positive (>=0)'))
	
	# units
	L = _globals.L_units[_geta(xobj, 'L. unit').string]
	F = _globals.F_units[_geta(xobj, 'F. unit').string]
	mm = L # to mm, from mm = 1/mm
	MPa = F/L/L # to MPa, from MPa = 1/MPa
	
	# convert inputs to N-m
	fc = fc*MPa
	cclear = cclear*mm
	
	# define basic points according to MC2020
	X, Y1, Y2, D, stress_null = _mc2020.make(pull_out, good_bond, confined, fc, cclear)
	
	# convert back to original units
	X = [ix/mm for ix in X]
	Y1 = [iy/MPa for iy in Y1]
	Y2 = [iy/MPa for iy in Y2]
	stress_null /= MPa
	
	# define null backbone points
	E1 = Y1[1]/X[1]
	E2 = Y2[1]/X[1]
	strain_null = stress_null/E1
	X_null = [0.0, strain_null, strain_null*2]
	Y_null = [0.0, stress_null, stress_null]
	D_null = [1.0, 1.0, 1.0]
	def to_tcl(x):
		return ' '.join(str(i) for i in x)
	
	# define 3 extra material IDs
	def next_id():
		id = pinfo.next_physicalProperties_id
		pinfo.next_physicalProperties_id += 1
		return id
	id_pos = next_id()
	id_neg = next_id()
	id_res = next_id()
	
	# write comment
	pinfo.out_file.write('{}# ASDBondSlip [{}]\n'.format(pinfo.indent, tag))
	
	# positive pinching material
	pinfo.out_file.write(_globals.concrete_format.format(
		pinfo.indent, id_pos, E1, eta, 
		to_tcl(X), to_tcl(Y1), to_tcl(D),
		to_tcl(X_null), to_tcl(Y_null), to_tcl(D_null),
		' -tangent' if tangent else '',
		' -implex' if implex else ''))
	# negative pinching material
	pinfo.out_file.write(_globals.concrete_format.format(
		pinfo.indent, id_neg, E1, eta, 
		to_tcl(X_null), to_tcl(Y_null), to_tcl(D_null),
		to_tcl(X), to_tcl(Y1), to_tcl(D),
		' -tangent' if tangent else '',
		' -implex' if implex else ''))
	# frictional material
	DP = [0.0]*len(X)
	pinfo.out_file.write(_globals.concrete_format.format(
		pinfo.indent, id_res, E2, eta, 
		to_tcl(X), to_tcl(Y2), to_tcl(DP),
		to_tcl(X), to_tcl(Y2), to_tcl(DP),
		' -tangent' if tangent else '',
		' -implex' if implex else ''))
	# combine
	pinfo.out_file.write('{}uniaxialMaterial Parallel {}  {} {} {} -factors {} {} {}\n'.format(pinfo.indent, tag, id_pos, id_neg, id_res, 1.0, 1.0, 1.0))
