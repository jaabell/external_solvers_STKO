import PyMpc.Units as u
from PyMpc import *
from mpc_utils_html import *
import PyMpc
import PyMpc.Math
import math
import opensees.utils.tcl_input as tclin

from scipy.signal import chirp, spectrogram
from scipy.fft import fft, ifft
from scipy.interpolate import CubicSpline
import numpy as np
#import eqsig

def makeXObjectMetaData():
    # Function
    at_Function = MpcAttributeMetaData()
    at_Function.type = MpcAttributeType.String
    at_Function.name = '__mpc_function__'
    at_Function.group = 'Group'
    at_Function.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Function')+'<br/>') + 
        html_par('') +
        html_par(html_href('','')+'<br/>') +
        html_end()
    )
    at_Function.editable = False

    # Select type of custom time series
    at_selectType = MpcAttributeMetaData()
    at_selectType.type = MpcAttributeType.String
    at_selectType.name = 'Select_Type'
    at_selectType.group = 'Group'
    at_selectType.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Select Type of Custom time series')+'<br/>') +
        html_par('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_selectType.sourceType = MpcAttributeSourceType.List
    at_selectType.setSourceList(['SineSweep', 'WhiteNoise', 'CubicSpline', 'Steps'])
    at_selectType.setDefault('SineSweep')

    # SineSweep
    at_sineSweep = MpcAttributeMetaData()
    at_sineSweep.type = MpcAttributeType.Boolean
    at_sineSweep.name = 'SineSweep'
    at_sineSweep.group = 'Group'
    at_sineSweep.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Sine Sweep')+'<br/>') +
        html_par('') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXX','Sine Sweep')+'<br/>') +
        html_end()
    )
    at_sineSweep.editable = False

    # f0
    at_f0 = MpcAttributeMetaData()
    at_f0.type = MpcAttributeType.Real
    at_f0.name = 'f0'
    at_f0.group = '-frequency'
    at_f0.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('f0')+'<br/>') +
        html_par('Initial frequency.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_f0.setDefault(2.0)

    # f1
    at_f1 = MpcAttributeMetaData()
    at_f1.type = MpcAttributeType.Real
    at_f1.name = 'f1'
    at_f1.group = '-frequency'
    at_f1.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('f1')+'<br/>') +
        html_par('Final frequency.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_f1.setDefault(1.0)

    # amplitude
    at_amplitude = MpcAttributeMetaData()
    at_amplitude.type = MpcAttributeType.Real
    at_amplitude.name = 'Amplitude'
    at_amplitude.group = '-Amplitude'
    at_amplitude.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Amplitude signal')+'<br/>') +
        html_par('Amplitude signal.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_amplitude.setDefault(0.1)

    # time of octaves
    at_timeForOctaves = MpcAttributeMetaData()
    at_timeForOctaves.type = MpcAttributeType.Real
    at_timeForOctaves.name = 'timeForOctaves'
    at_timeForOctaves.group = '-timeForOctaves'
    at_timeForOctaves.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Time of octaves')+'<br/>') +
        html_par('Time of octaves.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_timeForOctaves.setDefault(30.0)

    # Division
    at_Division = MpcAttributeMetaData()
    at_Division.type = MpcAttributeType.Real
    at_Division.name = 'Division'
    at_Division.group = '-Division'
    at_Division.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Division')+'<br/>') +
        html_par('Division.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_Division.setDefault(1000.0)

    # WhiteNoise
    at_WhiteNoise = MpcAttributeMetaData()
    at_WhiteNoise.type = MpcAttributeType.Boolean
    at_WhiteNoise.name = 'WhiteNoise'
    at_WhiteNoise.group = 'Group'
    at_WhiteNoise.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('White Noise')+'<br/>') +
        html_par('WhiteNoise') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_WhiteNoise.editable = False

    # ScaleF
    at_ScaleF = MpcAttributeMetaData()
    at_ScaleF.type = MpcAttributeType.Real
    at_ScaleF.name = 'ScaleF'
    at_ScaleF.group = '-Scale'
    at_ScaleF.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('ScaleF')+'<br/>') +
        html_par('Scale Factor.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_ScaleF.setDefault(1.0)

    # fmin
    at_fmin = MpcAttributeMetaData()
    at_fmin.type = MpcAttributeType.Real
    at_fmin.name = 'fmin'
    at_fmin.group = '-frequency'
    at_fmin.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('fmin')+'<br/>') +
        html_par('Minimum frequency.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_fmin.setDefault(100.0)

    # fmax
    at_fmax = MpcAttributeMetaData()
    at_fmax.type = MpcAttributeType.Real
    at_fmax.name = 'fmax'
    at_fmax.group = '-frequency'
    at_fmax.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('fmax')+'<br/>') +
        html_par('Maximum frequency.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_fmax.setDefault(500.0)

    # Sample Rate
    at_SampleRate = MpcAttributeMetaData()
    at_SampleRate.type = MpcAttributeType.Real
    at_SampleRate.name = 'SampleRate'
    at_SampleRate.group = '-Sample'
    at_SampleRate.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Sample Rate')+'<br/>') +
        html_par('Sample Rate.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_SampleRate.setDefault(2000.0)

    # Number of Samples
    at_NSamples = MpcAttributeMetaData()
    at_NSamples.type = MpcAttributeType.Integer
    at_NSamples.name = 'NSamples'
    at_NSamples.group = '-Sample'
    at_NSamples.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Number of Samples')+'<br/>') +
        html_par('Number of Samples.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_NSamples.setDefault(3500)

    # ===== CUBIC SPLINE ATTRIBUTES =====
    
    # CubicSpline
    at_CubicSpline = MpcAttributeMetaData()
    at_CubicSpline.type = MpcAttributeType.Boolean
    at_CubicSpline.name = 'CubicSpline'
    at_CubicSpline.group = 'Group'
    at_CubicSpline.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Cubic Spline')+'<br/>') +
        html_par('Cubic spline interpolation with configurable boundary conditions') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_CubicSpline.editable = False

    # Initial Time
    at_InitialTime = MpcAttributeMetaData()
    at_InitialTime.type = MpcAttributeType.Real
    at_InitialTime.name = 'InitialTime'
    at_InitialTime.group = '-Spline Parameters'
    at_InitialTime.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Initial Time')+'<br/>') +
        html_par('Start time for the cubic spline.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_InitialTime.setDefault(0.0)

    # Initial Value
    at_InitialValue = MpcAttributeMetaData()
    at_InitialValue.type = MpcAttributeType.Real
    at_InitialValue.name = 'InitialValue'
    at_InitialValue.group = '-Spline Parameters'
    at_InitialValue.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Initial Value')+'<br/>') +
        html_par('Value at the start time of the time series.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_InitialValue.setDefault(0.0)

    # Initial Slope
    at_InitialSlope = MpcAttributeMetaData()
    at_InitialSlope.type = MpcAttributeType.Real
    at_InitialSlope.name = 'InitialSlope'
    at_InitialSlope.group = '-Spline Parameters'
    at_InitialSlope.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Initial Slope')+'<br/>') +
        html_par('Derivative (slope) at the start of the time series.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_InitialSlope.setDefault(1.0)

    # Final Time
    at_FinalTime = MpcAttributeMetaData()
    at_FinalTime.type = MpcAttributeType.Real
    at_FinalTime.name = 'FinalTime'
    at_FinalTime.group = '-Spline Parameters'
    at_FinalTime.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Final Time')+'<br/>') +
        html_par('End time for the cubic spline.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_FinalTime.setDefault(10.0)

    # Final Value
    at_FinalValue = MpcAttributeMetaData()
    at_FinalValue.type = MpcAttributeType.Real
    at_FinalValue.name = 'FinalValue'
    at_FinalValue.group = '-Spline Parameters'
    at_FinalValue.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Final Value')+'<br/>') +
        html_par('Value at the end of the time series.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_FinalValue.setDefault(1.0)

    # Final Slope
    at_FinalSlope = MpcAttributeMetaData()
    at_FinalSlope.type = MpcAttributeType.Real
    at_FinalSlope.name = 'FinalSlope'
    at_FinalSlope.group = '-Spline Parameters'
    at_FinalSlope.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Final Slope')+'<br/>') +
        html_par('Derivative (slope) at the end of the time series.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_FinalSlope.setDefault(0.0)

    # Number of Points for Spline
    at_SplinePoints = MpcAttributeMetaData()
    at_SplinePoints.type = MpcAttributeType.Integer
    at_SplinePoints.name = 'SplinePoints'
    at_SplinePoints.group = '-Spline Parameters'
    at_SplinePoints.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Number of Points')+'<br/>') +
        html_par('Number of points for spline evaluation.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_SplinePoints.setDefault(1000)

    # ===== STEPS ATTRIBUTES =====
    
    # Steps
    at_Steps = MpcAttributeMetaData()
    at_Steps.type = MpcAttributeType.Boolean
    at_Steps.name = 'Steps'
    at_Steps.group = 'Group'
    at_Steps.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Steps')+'<br/>') +
        html_par('Step function with configurable step length and height') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_Steps.editable = False

    # Start Time
    at_StartTime = MpcAttributeMetaData()
    at_StartTime.type = MpcAttributeType.Real
    at_StartTime.name = 'StartTime'
    at_StartTime.group = '-Step Parameters'
    at_StartTime.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Start Time')+'<br/>') +
        html_par('Start time for the step function.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_StartTime.setDefault(0.0)

    # End Time
    at_EndTime = MpcAttributeMetaData()
    at_EndTime.type = MpcAttributeType.Real
    at_EndTime.name = 'EndTime'
    at_EndTime.group = '-Step Parameters'
    at_EndTime.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('End Time')+'<br/>') +
        html_par('End time for the step function.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_EndTime.setDefault(10.0)

    # Step Length
    at_StepLength = MpcAttributeMetaData()
    at_StepLength.type = MpcAttributeType.Real
    at_StepLength.name = 'StepLength'
    at_StepLength.group = '-Step Parameters'
    at_StepLength.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Step Length')+'<br/>') +
        html_par('Duration of each step in time units.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_StepLength.setDefault(1.0)

    # Step Height
    at_StepHeight = MpcAttributeMetaData()
    at_StepHeight.type = MpcAttributeType.Real
    at_StepHeight.name = 'StepHeight'
    at_StepHeight.group = '-Step Parameters'
    at_StepHeight.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Step Height')+'<br/>') +
        html_par('Height increase for each step.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_StepHeight.setDefault(1.0)

    # Transition Interval
    at_TransitionInterval = MpcAttributeMetaData()
    at_TransitionInterval.type = MpcAttributeType.Real
    at_TransitionInterval.name = 'TransitionInterval'
    at_TransitionInterval.group = '-Step Parameters'
    at_TransitionInterval.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Transition Interval')+'<br/>') +
        html_par('Small time interval for the step transition to maintain monotonic time vector.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_TransitionInterval.setDefault(1e-6)

    # Initial Value for Steps
    at_StepInitialValue = MpcAttributeMetaData()
    at_StepInitialValue.type = MpcAttributeType.Real
    at_StepInitialValue.name = 'StepInitialValue'
    at_StepInitialValue.group = '-Step Parameters'
    at_StepInitialValue.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Initial Value')+'<br/>') +
        html_par('Initial value before the first step.') +
        html_par(html_href('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX','Custom TimeSeries')+'<br/>') +
        html_end()
    )
    at_StepInitialValue.setDefault(0.0)

    # PLOT OPTION
    # Data type
    at_DataType = MpcAttributeMetaData()
    at_DataType.type = MpcAttributeType.String
    at_DataType.name = 'Data_Type_Plot'
    at_DataType.group = 'Plot Options'
    at_DataType.description = (
        html_par(html_begin()) +
        html_par(html_boldtext('Data type')+'<br/>') +
        html_par('Data type') +
        html_end()
    )
    at_DataType.sourceType = MpcAttributeSourceType.List
    at_DataType.setSourceList(['TimeSeries', 'FourierSpectra', 'ResponceSpectra'])
    at_DataType.setDefault('TimeSeries')

    xom = MpcXObjectMetaData()
    xom.name = 'CustomTimeSeries'
    xom.addAttribute(at_Function)
    xom.addAttribute(at_selectType)
    xom.addAttribute(at_sineSweep)
    xom.addAttribute(at_f0)
    xom.addAttribute(at_f1)
    xom.addAttribute(at_amplitude)
    xom.addAttribute(at_timeForOctaves)
    xom.addAttribute(at_Division)
    xom.addAttribute(at_WhiteNoise)
    xom.addAttribute(at_ScaleF)
    xom.addAttribute(at_fmin)
    xom.addAttribute(at_fmax)
    xom.addAttribute(at_SampleRate)
    xom.addAttribute(at_NSamples)
    xom.addAttribute(at_CubicSpline)
    xom.addAttribute(at_InitialTime)
    xom.addAttribute(at_InitialValue)
    xom.addAttribute(at_InitialSlope)
    xom.addAttribute(at_FinalTime)
    xom.addAttribute(at_FinalValue)
    xom.addAttribute(at_FinalSlope)
    xom.addAttribute(at_SplinePoints)
    xom.addAttribute(at_Steps)
    xom.addAttribute(at_StartTime)
    xom.addAttribute(at_EndTime)
    xom.addAttribute(at_StepLength)
    xom.addAttribute(at_StepHeight)
    xom.addAttribute(at_TransitionInterval)
    xom.addAttribute(at_StepInitialValue)
    xom.addAttribute(at_DataType)

    # Sine Sweep dependencies
    xom.setVisibilityDependency(at_sineSweep, at_f0)
    xom.setVisibilityDependency(at_sineSweep, at_f1)
    xom.setVisibilityDependency(at_sineSweep, at_amplitude)
    xom.setVisibilityDependency(at_sineSweep, at_timeForOctaves)
    xom.setVisibilityDependency(at_sineSweep, at_Division)

    # White Noise dependencies
    xom.setVisibilityDependency(at_WhiteNoise, at_fmin)
    xom.setVisibilityDependency(at_WhiteNoise, at_fmax)
    xom.setVisibilityDependency(at_WhiteNoise, at_SampleRate)
    xom.setVisibilityDependency(at_WhiteNoise, at_NSamples)
    xom.setVisibilityDependency(at_WhiteNoise, at_ScaleF)

    # Cubic Spline dependencies
    xom.setVisibilityDependency(at_CubicSpline, at_InitialTime)
    xom.setVisibilityDependency(at_CubicSpline, at_InitialValue)
    xom.setVisibilityDependency(at_CubicSpline, at_InitialSlope)
    xom.setVisibilityDependency(at_CubicSpline, at_FinalTime)
    xom.setVisibilityDependency(at_CubicSpline, at_FinalValue)
    xom.setVisibilityDependency(at_CubicSpline, at_FinalSlope)
    xom.setVisibilityDependency(at_CubicSpline, at_SplinePoints)

    # Steps dependencies
    xom.setVisibilityDependency(at_Steps, at_StartTime)
    xom.setVisibilityDependency(at_Steps, at_EndTime)
    xom.setVisibilityDependency(at_Steps, at_StepLength)
    xom.setVisibilityDependency(at_Steps, at_StepHeight)
    xom.setVisibilityDependency(at_Steps, at_TransitionInterval)
    xom.setVisibilityDependency(at_Steps, at_StepInitialValue)

    # auto-exclusive dependencies
    xom.setBooleanAutoExclusiveDependency(at_selectType, at_sineSweep)
    xom.setBooleanAutoExclusiveDependency(at_selectType, at_WhiteNoise)
    xom.setBooleanAutoExclusiveDependency(at_selectType, at_CubicSpline)
    xom.setBooleanAutoExclusiveDependency(at_selectType, at_Steps)

    return xom

def evaluateFunctionAttribute(xobj):
    if(xobj is None):
        print('Error: xobj is null\n')
        return PyMpc.Math.mat()

    if(xobj.name != 'CustomTimeSeries'):
        print('Error: invalid xobj type, expected "CustomTimeSeries", given "{}"'.format(xobj.name))
        return PyMpc.Math.mat()

    func_at = xobj.getAttribute('__mpc_function__')
    if(func_at is None):
        print('Error: cannot find "__mpc_function__" attribute\n')
        return PyMpc.Math.mat()

    selectType_at = xobj.getAttribute('Select_Type')
    if(selectType_at is None):
        print('Error: cannot find "Select_Type" attribute\n')
        return PyMpc.Math.mat()

    SineSweep_at = xobj.getAttribute('SineSweep')
    if(SineSweep_at is None):
        print('Error: cannot find "SineSweep" attribute\n')
        return PyMpc.Math.mat()

    WhiteNoise_at = xobj.getAttribute('WhiteNoise')
    if(WhiteNoise_at is None):
        print('Error: cannot find "WhiteNoise" attribute\n')
        return PyMpc.Math.mat()

    CubicSpline_at = xobj.getAttribute('CubicSpline')
    if(CubicSpline_at is None):
        print('Error: cannot find "CubicSpline" attribute\n')
        return PyMpc.Math.mat()

    Steps_at = xobj.getAttribute('Steps')
    if(Steps_at is None):
        print('Error: cannot find "Steps" attribute\n')
        return PyMpc.Math.mat()

    DataTypePlot_at = xobj.getAttribute('Data_Type_Plot')
    if(DataTypePlot_at is None):
        print('Error: cannot find "Data_Type_Plot" attribute\n')
        return PyMpc.Math.mat()

    DataTypePlot_at = DataTypePlot_at.string
    SineSweep = SineSweep_at.boolean
    WhiteNoise = WhiteNoise_at.boolean
    UseCubicSpline = CubicSpline_at.boolean
    Steps_at = xobj.getAttribute('Steps')
    if(Steps_at is None):
        raise Exception('Error: cannot find "Steps" attribute')
    UseSteps = Steps_at.boolean

    if(SineSweep == True):
        f0_at = xobj.getAttribute('f0')
        if(f0_at is None):
            print('Error: cannot find "f0" attribute\n')
            return PyMpc.Math.mat()
        f1_at = xobj.getAttribute('f1')
        if(f1_at is None):
            print('Error: cannot find "f1" attribute\n')
            return PyMpc.Math.mat()
        Amplitude_at = xobj.getAttribute('Amplitude')
        if(Amplitude_at is None):
            print('Error: cannot find "Amplitude" attribute\n')
            return PyMpc.Math.mat()
        timeForOctaves_at = xobj.getAttribute('timeForOctaves')
        if(timeForOctaves_at is None):
            print('Error: cannot find "timeForOctaves" attribute\n')
            return PyMpc.Math.mat()
        Division_at = xobj.getAttribute('Division')
        if(Division_at is None):
            print('Error: cannot find "Division" attribute\n')
            return PyMpc.Math.mat()

        f0=f0_at.real
        f1=f1_at.real
        Amplitude=Amplitude_at.real
        timeForOctaves=timeForOctaves_at.real
        Division=Division_at.real

        NOctave=int(np.log(f0/f1)/np.log(2))
        t1=timeForOctaves*NOctave
        f0_last=f1
        for i in range(0,NOctave):
            f0_last=f0_last*2
        f0=f0_last

        t = np.linspace(0, timeForOctaves*NOctave, int(timeForOctaves*NOctave*Division))
        x = chirp(t, f0=f0, f1=f1, t1=t1, method='linear')*Amplitude

        xy = PyMpc.Math.mat(len(x), 2)
        for i in range(len(t)):
            xy[i, 0] = t[i]
            xy[i, 1] = x[i]

    elif(WhiteNoise == True):
        fmin_at = xobj.getAttribute('fmin')
        if(fmin_at is None):
            print('Error: cannot find "fmin" attribute\n')
            return PyMpc.Math.mat()
        fmax_at = xobj.getAttribute('fmax')
        if(fmax_at is None):
            print('Error: cannot find "fmax" attribute\n')
            return PyMpc.Math.mat()
        SampleRate_at = xobj.getAttribute('SampleRate')
        if(SampleRate_at is None):
            print('Error: cannot find "SampleRate" attribute\n')
            return PyMpc.Math.mat()
        NSamples_at = xobj.getAttribute('NSamples')
        if(NSamples_at is None):
            print('Error: cannot find "NSamples" attribute\n')
            return PyMpc.Math.mat()
        ScaleF_at = xobj.getAttribute('ScaleF')
        if(ScaleF_at is None):
            print('Error: cannot find "ScaleF" attribute\n')
            return PyMpc.Math.mat()

        def fftnoise(f):
            f = np.array(f, dtype='complex')
            Np = (len(f) - 1) // 2
            phases = np.random.rand(Np) * 2 * np.pi
            phases = np.cos(phases) + 1j * np.sin(phases)
            f[1:Np+1] *= phases
            f[-1:-1-Np:-1] = np.conj(f[1:Np+1])
            return np.fft.ifft(f).real

        def band_limited_noise(min_freq, max_freq, samples, samplerate):
            freqs = np.abs(np.fft.fftfreq(samples, 1/samplerate))
            f = np.zeros(samples)
            idx = np.where(np.logical_and(freqs>=min_freq, freqs<=max_freq))[0]
            f[idx] = 1
            return fftnoise(f)

        ScaleF=ScaleF_at.real
        Fmin=fmin_at.real
        Fmax=fmax_at.real
        SampleRate=SampleRate_at.real
        Sample=NSamples_at.integer

        Durate=Sample*(1/SampleRate)
        t=np.linspace(0.0,Durate,int(SampleRate*(Sample/SampleRate)))
        x = band_limited_noise(Fmin, Fmax, Sample, SampleRate)
        x = x * (2**15 - 1)

        xy = PyMpc.Math.mat(len(x), 2)
        x=list(x)
        for i in range(len(x)):
            xy[i, 0] = t[i]
            xy[i, 1] = x[i]*ScaleF

    elif(UseCubicSpline == True):
        # Get cubic spline parameters
        InitialTime_at = xobj.getAttribute('InitialTime')
        if(InitialTime_at is None):
            print('Error: cannot find "InitialTime" attribute\n')
            return PyMpc.Math.mat()
        InitialValue_at = xobj.getAttribute('InitialValue')
        if(InitialValue_at is None):
            print('Error: cannot find "InitialValue" attribute\n')
            return PyMpc.Math.mat()
        InitialSlope_at = xobj.getAttribute('InitialSlope')
        if(InitialSlope_at is None):
            print('Error: cannot find "InitialSlope" attribute\n')
            return PyMpc.Math.mat()
        FinalTime_at = xobj.getAttribute('FinalTime')
        if(FinalTime_at is None):
            print('Error: cannot find "FinalTime" attribute\n')
            return PyMpc.Math.mat()
        FinalValue_at = xobj.getAttribute('FinalValue')
        if(FinalValue_at is None):
            print('Error: cannot find "FinalValue" attribute\n')
            return PyMpc.Math.mat()
        FinalSlope_at = xobj.getAttribute('FinalSlope')
        if(FinalSlope_at is None):
            print('Error: cannot find "FinalSlope" attribute\n')
            return PyMpc.Math.mat()
        SplinePoints_at = xobj.getAttribute('SplinePoints')
        if(SplinePoints_at is None):
            print('Error: cannot find "SplinePoints" attribute\n')
            return PyMpc.Math.mat()

        # Extract values
        initial_time = InitialTime_at.real
        initial_value = InitialValue_at.real
        initial_slope = InitialSlope_at.real
        final_time = FinalTime_at.real
        final_value = FinalValue_at.real
        final_slope = FinalSlope_at.real
        n_points = SplinePoints_at.integer

        # Validate inputs
        if final_time <= initial_time:
            print('Error: Final time must be greater than initial time\n')
            return PyMpc.Math.mat()
        if n_points < 2:
            print('Error: Number of points must be at least 2\n')
            return PyMpc.Math.mat()

        # Create time vector
        t = np.linspace(initial_time, final_time, n_points)
        
        # Define boundary conditions for cubic spline
        # We need at least 2 points to define a spline with boundary conditions
        # Let's use the initial and final points with their slopes
        x_key = np.array([initial_time, final_time])
        y_key = np.array([initial_value, final_value])
        
        # Create cubic spline with clamped boundary conditions (specified derivatives)
        cs = CubicSpline(x_key, y_key, bc_type=((1, initial_slope), (1, final_slope)))
        
        # Evaluate the spline
        x = cs(t)

        # Create the output matrix
        xy = PyMpc.Math.mat(len(t), 2)
        for i in range(len(t)):
            xy[i, 0] = t[i]
            xy[i, 1] = x[i]

    elif(UseSteps == True):
        # Get steps parameters
        StartTime_at = xobj.getAttribute('StartTime')
        if(StartTime_at is None):
            print('Error: cannot find "StartTime" attribute\n')
            return PyMpc.Math.mat()
        EndTime_at = xobj.getAttribute('EndTime')
        if(EndTime_at is None):
            print('Error: cannot find "EndTime" attribute\n')
            return PyMpc.Math.mat()
        StepLength_at = xobj.getAttribute('StepLength')
        if(StepLength_at is None):
            print('Error: cannot find "StepLength" attribute\n')
            return PyMpc.Math.mat()
        StepHeight_at = xobj.getAttribute('StepHeight')
        if(StepHeight_at is None):
            print('Error: cannot find "StepHeight" attribute\n')
            return PyMpc.Math.mat()
        TransitionInterval_at = xobj.getAttribute('TransitionInterval')
        if(TransitionInterval_at is None):
            print('Error: cannot find "TransitionInterval" attribute\n')
            return PyMpc.Math.mat()
        StepInitialValue_at = xobj.getAttribute('StepInitialValue')
        if(StepInitialValue_at is None):
            print('Error: cannot find "StepInitialValue" attribute\n')
            return PyMpc.Math.mat()

        # Extract values
        start_time = StartTime_at.real
        end_time = EndTime_at.real
        step_length = StepLength_at.real
        step_height = StepHeight_at.real
        transition_interval = TransitionInterval_at.real
        initial_value = StepInitialValue_at.real

        # Validate inputs
        if end_time <= start_time:
            print('Error: End time must be greater than start time\n')
            return PyMpc.Math.mat()
        if step_length <= 0:
            print('Error: Step length must be positive\n')
            return PyMpc.Math.mat()
        if transition_interval <= 0:
            print('Error: Transition interval must be positive\n')
            return PyMpc.Math.mat()
        if transition_interval >= step_length:
            print('Error: Transition interval must be smaller than step length\n')
            return PyMpc.Math.mat()

        # Generate step function
        total_duration = end_time - start_time
        num_steps = int(np.floor(total_duration / step_length))
        
        # Create time and value arrays
        t_list = []
        x_list = []
        
        # Add initial point
        t_list.append(start_time)
        x_list.append(initial_value)
        
        current_value = initial_value
        current_time = start_time
        
        for step in range(num_steps):
            # Time just before the step
            step_time = start_time + step * step_length
            if step_time > current_time:
                t_list.append(step_time)
                x_list.append(current_value)
            
            # Time just after the step (with transition interval)
            step_time_after = step_time + transition_interval
            current_value += step_height
            t_list.append(step_time_after)
            x_list.append(current_value)
            
            current_time = step_time_after
        
        # Add final point if needed
        if current_time < end_time:
            t_list.append(end_time)
            x_list.append(current_value)
        
        # Convert to numpy arrays
        t = np.array(t_list)
        x = np.array(x_list)

        # Create the output matrix
        xy = PyMpc.Math.mat(len(t), 2)
        for i in range(len(t)):
            xy[i, 0] = t[i]
            xy[i, 1] = x[i]

    # if DataTypePlot_at=='FourierSpectra':
    #     org_signal = eqsig.Signal(x, t[1]-t[0])
    #     Tmax=5.0
    #     d_step = Tmax / 100
    #     T = np.arange(d_step * 1, Tmax, d_step)
    #     org_signal.gen_fa_spectrum()
    #     FuSpec=abs(org_signal.fa_spectrum)
    #     FuFreq=org_signal.fa_frequencies
    #     xy = PyMpc.Math.mat(len(FuFreq), 2)
    #     for i in range(len(FuFreq)):
    #         xy[i, 0] = FuFreq[i]
    #         xy[i, 1] = FuSpec[i]

    # if DataTypePlot_at=='ResponceSpectra':
    #     org_signal = eqsig.AccSignal(x, t[1]-t[0])
    #     Tmax=5.0
    #     d_step = Tmax / 100
    #     T = np.arange(d_step * 1, Tmax, d_step)
    #     xy = PyMpc.Math.mat(len(T), 2)
    #     org_signal.generate_response_spectrum(response_times=T,xi=0.05)
    #     SpSc=org_signal.s_a
    #     for i in range(len(SpSc)):
    #         xy[i, 0] = T[i]
    #         xy[i, 1] = SpSc[i]

    return xy

def writeTcl(pinfo):
    xobj = pinfo.definition.XObject
    ClassName = xobj.name
    if pinfo.currentDescription != ClassName:
        pinfo.out_file.write('\n{}# {} {}\n'.format(pinfo.indent, xobj.Xnamespace, ClassName))
        pinfo.currentDescription = ClassName

    tag = xobj.parent.componentId
    sopt = ''  # optional string

    # Get common attributes
    SineSweep_at = xobj.getAttribute('SineSweep')
    if(SineSweep_at is None):
        raise Exception('Error: cannot find "SineSweep" attribute')
    SineSweep = SineSweep_at.boolean

    WhiteNoise_at = xobj.getAttribute('WhiteNoise')
    if(WhiteNoise_at is None):
        raise Exception('Error: cannot find "WhiteNoise" attribute')
    WhiteNoise = WhiteNoise_at.boolean

    CubicSpline_at = xobj.getAttribute('CubicSpline')
    if(CubicSpline_at is None):
        raise Exception('Error: cannot find "CubicSpline" attribute')
    UseCubicSpline = CubicSpline_at.boolean

    Steps_at = xobj.getAttribute('Steps')
    if(Steps_at is None):
        raise Exception('Error: cannot find "Steps" attribute')
    UseSteps = Steps_at.boolean

    if SineSweep:
        # SineSweep implementation
        f0_at = xobj.getAttribute('f0')
        if(f0_at is None):
            raise Exception('Error: cannot find "f0" attribute')
        f0 = f0_at.real

        f1_at = xobj.getAttribute('f1')
        if(f1_at is None):
            raise Exception('Error: cannot find "f1" attribute')
        f1 = f1_at.real

        Amplitude_at = xobj.getAttribute('Amplitude')
        if(Amplitude_at is None):
            raise Exception('Error: cannot find "Amplitude" attribute')
        Amplitude = Amplitude_at.real

        timeForOctaves_at = xobj.getAttribute('timeForOctaves')
        if(timeForOctaves_at is None):
            raise Exception('Error: cannot find "timeForOctaves" attribute')
        timeForOctaves = timeForOctaves_at.real

        Division_at = xobj.getAttribute('Division')
        if(Division_at is None):
            raise Exception('Error: cannot find "Division" attribute')
        Division = Division_at.real

        NOctave=int(np.log(f0/f1)/np.log(2))
        t1=timeForOctaves*NOctave
        f0_last=f1
        for i in range(0,NOctave):
            f0_last=f0_last*2
        f0=f0_last

        t = np.linspace(0, timeForOctaves*NOctave, int(timeForOctaves*NOctave*Division))
        w = chirp(t, f0=f0, f1=f1, t1=t1, method='linear')*Amplitude

        listTimes = t
        listValues = w

        # Write TCL lists
        times_str = 'set timeSeries_list_of_times_{}'.format(tag)+' {'
        values_str = 'set timeSeries_list_of_values_{}'.format(tag)+' {'
        nLettersT = len(times_str)
        nLettersV = len(values_str)
        nTabT = nLettersT // 4
        nTabV = nLettersV // 4
        n = 1
        for i in range(len(listTimes)):
            if (i == (10*n)):
                times_str += '\\\n{}{}'.format(pinfo.indent, tclin.utils.nIndent(nTabT))
                values_str += '\\\n{}{}'.format(pinfo.indent, tclin.utils.nIndent(nTabV))
                n += 1
            if (i!=len(listTimes)-1):
                times_str += '{} '.format(listTimes[i])
                values_str += '{} '.format(listValues[i])
            else:
                times_str += '{}'.format(listTimes[i])
                values_str += '{}'.format(listValues[i])
        times_str += '}\n'
        values_str += '}\n'

        pinfo.out_file.write(times_str)
        pinfo.out_file.write(values_str)

        str_tcl = '{0}timeSeries Path {1} -time $timeSeries_list_of_times_{1} -values $timeSeries_list_of_values_{1}{2}\n'.format(pinfo.indent, tag, sopt)

    elif UseSteps:
        # Steps implementation
        StartTime_at = xobj.getAttribute('StartTime')
        if(StartTime_at is None):
            raise Exception('Error: cannot find "StartTime" attribute')
        start_time = StartTime_at.real

        EndTime_at = xobj.getAttribute('EndTime')
        if(EndTime_at is None):
            raise Exception('Error: cannot find "EndTime" attribute')
        end_time = EndTime_at.real

        StepLength_at = xobj.getAttribute('StepLength')
        if(StepLength_at is None):
            raise Exception('Error: cannot find "StepLength" attribute')
        step_length = StepLength_at.real

        StepHeight_at = xobj.getAttribute('StepHeight')
        if(StepHeight_at is None):
            raise Exception('Error: cannot find "StepHeight" attribute')
        step_height = StepHeight_at.real

        TransitionInterval_at = xobj.getAttribute('TransitionInterval')
        if(TransitionInterval_at is None):
            raise Exception('Error: cannot find "TransitionInterval" attribute')
        transition_interval = TransitionInterval_at.real

        StepInitialValue_at = xobj.getAttribute('StepInitialValue')
        if(StepInitialValue_at is None):
            raise Exception('Error: cannot find "StepInitialValue" attribute')
        initial_value = StepInitialValue_at.real

        # Validate inputs
        if end_time <= start_time:
            raise Exception('Error: End time must be greater than start time')
        if step_length <= 0:
            raise Exception('Error: Step length must be positive')
        if transition_interval <= 0:
            raise Exception('Error: Transition interval must be positive')
        if transition_interval >= step_length:
            raise Exception('Error: Transition interval must be smaller than step length')

        # Generate step function
        total_duration = end_time - start_time
        num_steps = int(np.floor(total_duration / step_length))
        
        # Create time and value arrays
        t_list = []
        x_list = []
        
        # Add initial point
        t_list.append(start_time)
        x_list.append(initial_value)
        
        current_value = initial_value
        current_time = start_time
        
        for step in range(num_steps):
            # Time just before the step
            step_time = start_time + step * step_length
            if step_time > current_time:
                t_list.append(step_time)
                x_list.append(current_value)
            
            # Time just after the step (with transition interval)
            step_time_after = step_time + transition_interval
            current_value += step_height
            t_list.append(step_time_after)
            x_list.append(current_value)
            
            current_time = step_time_after
        
        # Add final point if needed
        if current_time < end_time:
            t_list.append(end_time)
            x_list.append(current_value)
        
        # Convert to numpy arrays
        t = np.array(t_list)
        x = np.array(x_list)

        listTimes = t
        listValues = x

        # Write TCL lists
        times_str = 'set timeSeries_list_of_times_{}'.format(tag)+' {'
        values_str = 'set timeSeries_list_of_values_{}'.format(tag)+' {'
        nLettersT = len(times_str)
        nLettersV = len(values_str)
        nTabT = nLettersT // 4
        nTabV = nLettersV // 4
        n = 1
        for i in range(len(listTimes)):
            if (i == (10*n)):
                times_str += '\\\n{}{}'.format(pinfo.indent, tclin.utils.nIndent(nTabT))
                values_str += '\\\n{}{}'.format(pinfo.indent, tclin.utils.nIndent(nTabV))
                n += 1
            if (i!=len(listTimes)-1):
                times_str += '{} '.format(listTimes[i])
                values_str += '{} '.format(listValues[i])
            else:
                times_str += '{}'.format(listTimes[i])
                values_str += '{}'.format(listValues[i])
        times_str += '}\n'
        values_str += '}\n'

        pinfo.out_file.write(times_str)
        pinfo.out_file.write(values_str)

        str_tcl = '{0}timeSeries Path {1} -time $timeSeries_list_of_times_{1} -values $timeSeries_list_of_values_{1}{2}\n'.format(pinfo.indent, tag, sopt)

    elif WhiteNoise:
        # WhiteNoise implementation
        ScaleF_at = xobj.getAttribute('ScaleF')
        if(ScaleF_at is None):
            raise Exception('Error: cannot find "ScaleF" attribute')
        ScaleF = ScaleF_at.real

        fmin_at = xobj.getAttribute('fmin')
        if(fmin_at is None):
            raise Exception('Error: cannot find "fmin" attribute')
        Fmin = fmin_at.real

        fmax_at = xobj.getAttribute('fmax')
        if(fmax_at is None):
            raise Exception('Error: cannot find "fmax" attribute')
        Fmax = fmax_at.real

        SampleRate_at = xobj.getAttribute('SampleRate')
        if(SampleRate_at is None):
            raise Exception('Error: cannot find "SampleRate" attribute')
        SampleRate = SampleRate_at.real

        NSamples_at = xobj.getAttribute('NSamples')
        if(NSamples_at is None):
            raise Exception('Error: cannot find "NSamples" attribute')
        Sample = NSamples_at.integer

        def fftnoise(f):
            f = np.array(f, dtype='complex')
            Np = (len(f) - 1) // 2
            phases = np.random.rand(Np) * 2 * np.pi
            phases = np.cos(phases) + 1j * np.sin(phases)
            f[1:Np+1] *= phases
            f[-1:-1-Np:-1] = np.conj(f[1:Np+1])
            return np.fft.ifft(f).real

        def band_limited_noise(min_freq, max_freq, samples, samplerate):
            freqs = np.abs(np.fft.fftfreq(samples, 1/samplerate))
            f = np.zeros(samples)
            idx = np.where(np.logical_and(freqs>=min_freq, freqs<=max_freq))[0]
            f[idx] = 1
            return fftnoise(f)

        Durate=Sample*(1/SampleRate)
        t=np.linspace(0.0,Durate,int(SampleRate*(Sample/SampleRate)))
        x = band_limited_noise(Fmin, Fmax, Sample, SampleRate)
        x = x * (2**15 - 1)
        x=[xi*ScaleF for xi in x]

        listTimes = t
        listValues = x

        # Write TCL lists
        times_str = 'set timeSeries_list_of_times_{}'.format(tag)+' {'
        values_str = 'set timeSeries_list_of_values_{}'.format(tag)+' {'
        nLettersT = len(times_str)
        nLettersV = len(values_str)
        nTabT = nLettersT // 4
        nTabV = nLettersV // 4
        n = 1
        for i in range(len(listTimes)):
            if (i == (10*n)):
                times_str += '\\\n{}{}'.format(pinfo.indent, tclin.utils.nIndent(nTabT))
                values_str += '\\\n{}{}'.format(pinfo.indent, tclin.utils.nIndent(nTabV))
                n += 1
            if (i!=len(listTimes)-1):
                times_str += '{} '.format(listTimes[i])
                values_str += '{} '.format(listValues[i])
            else:
                times_str += '{}'.format(listTimes[i])
                values_str += '{}'.format(listValues[i])
        times_str += '}\n'
        values_str += '}\n'

        pinfo.out_file.write(times_str)
        pinfo.out_file.write(values_str)

        str_tcl = '{0}timeSeries Path {1} -time $timeSeries_list_of_times_{1} -values $timeSeries_list_of_values_{1}{2}\n'.format(pinfo.indent, tag, sopt)

    elif UseCubicSpline:
        # CubicSpline implementation
        InitialTime_at = xobj.getAttribute('InitialTime')
        if(InitialTime_at is None):
            raise Exception('Error: cannot find "InitialTime" attribute')
        initial_time = InitialTime_at.real

        InitialValue_at = xobj.getAttribute('InitialValue')
        if(InitialValue_at is None):
            raise Exception('Error: cannot find "InitialValue" attribute')
        initial_value = InitialValue_at.real

        InitialSlope_at = xobj.getAttribute('InitialSlope')
        if(InitialSlope_at is None):
            raise Exception('Error: cannot find "InitialSlope" attribute')
        initial_slope = InitialSlope_at.real

        FinalTime_at = xobj.getAttribute('FinalTime')
        if(FinalTime_at is None):
            raise Exception('Error: cannot find "FinalTime" attribute')
        final_time = FinalTime_at.real

        FinalValue_at = xobj.getAttribute('FinalValue')
        if(FinalValue_at is None):
            raise Exception('Error: cannot find "FinalValue" attribute')
        final_value = FinalValue_at.real

        FinalSlope_at = xobj.getAttribute('FinalSlope')
        if(FinalSlope_at is None):
            raise Exception('Error: cannot find "FinalSlope" attribute')
        final_slope = FinalSlope_at.real

        SplinePoints_at = xobj.getAttribute('SplinePoints')
        if(SplinePoints_at is None):
            raise Exception('Error: cannot find "SplinePoints" attribute')
        n_points = SplinePoints_at.integer

        # Validate inputs
        if final_time <= initial_time:
            raise Exception('Error: Final time must be greater than initial time')
        if n_points < 2:
            raise Exception('Error: Number of points must be at least 2')

        # Create time vector
        t = np.linspace(initial_time, final_time, n_points)
        
        # Define boundary conditions for cubic spline
        x_key = np.array([initial_time, final_time])
        y_key = np.array([initial_value, final_value])
        
        # Create cubic spline with clamped boundary conditions
        cs = CubicSpline(x_key, y_key, bc_type=((1, initial_slope), (1, final_slope)))
        
        # Evaluate the spline
        x = cs(t)

        listTimes = t
        listValues = x

        # Write TCL lists
        times_str = 'set timeSeries_list_of_times_{}'.format(tag)+' {'
        values_str = 'set timeSeries_list_of_values_{}'.format(tag)+' {'
        nLettersT = len(times_str)
        nLettersV = len(values_str)
        nTabT = nLettersT // 4
        nTabV = nLettersV // 4
        n = 1
        for i in range(len(listTimes)):
            if (i == (10*n)):
                times_str += '\\\n{}{}'.format(pinfo.indent, tclin.utils.nIndent(nTabT))
                values_str += '\\\n{}{}'.format(pinfo.indent, tclin.utils.nIndent(nTabV))
                n += 1
            if (i!=len(listTimes)-1):
                times_str += '{} '.format(listTimes[i])
                values_str += '{} '.format(listValues[i])
            else:
                times_str += '{}'.format(listTimes[i])
                values_str += '{}'.format(listValues[i])
        times_str += '}\n'
        values_str += '}\n'

        pinfo.out_file.write(times_str)
        pinfo.out_file.write(values_str)

        str_tcl = '{0}timeSeries Path {1} -time $timeSeries_list_of_times_{1} -values $timeSeries_list_of_values_{1}{2}\n'.format(pinfo.indent, tag, sopt)

    pinfo.out_file.write(str_tcl)