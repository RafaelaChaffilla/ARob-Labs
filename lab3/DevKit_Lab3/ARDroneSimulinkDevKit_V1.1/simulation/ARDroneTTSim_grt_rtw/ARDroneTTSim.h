/*
 * ARDroneTTSim.h
 *
 * Code generation for model "ARDroneTTSim".
 *
 * Model version              : $Id: UAV_SIL.mdl 912 2012-10-10 23:11:42Z joh07594 $
 * Simulink Coder version : 9.1 (R2019a) 23-Nov-2018
 * C source code generated on : Wed Dec 21 00:02:56 2022
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Freescale->32-bit PowerPC
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_ARDroneTTSim_h_
#define RTW_HEADER_ARDroneTTSim_h_
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>
#ifndef ARDroneTTSim_COMMON_INCLUDES_
# define ARDroneTTSim_COMMON_INCLUDES_
#include <stdlib.h>
#include "rtwtypes.h"
#include "zero_crossing_types.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#endif                                 /* ARDroneTTSim_COMMON_INCLUDES_ */

#include "ARDroneTTSim_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlockIO
# define rtmGetBlockIO(rtm)            ((rtm)->blockIO)
#endif

#ifndef rtmSetBlockIO
# define rtmSetBlockIO(rtm, val)       ((rtm)->blockIO = (val))
#endif

#ifndef rtmGetChecksums
# define rtmGetChecksums(rtm)          ((rtm)->Sizes.checksums)
#endif

#ifndef rtmSetChecksums
# define rtmSetChecksums(rtm, val)     ((rtm)->Sizes.checksums = (val))
#endif

#ifndef rtmGetConstBlockIO
# define rtmGetConstBlockIO(rtm)       ((rtm)->constBlockIO)
#endif

#ifndef rtmSetConstBlockIO
# define rtmSetConstBlockIO(rtm, val)  ((rtm)->constBlockIO = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetCtrlRateMdlRefTiming
# define rtmGetCtrlRateMdlRefTiming(rtm) ()
#endif

#ifndef rtmSetCtrlRateMdlRefTiming
# define rtmSetCtrlRateMdlRefTiming(rtm, val) ()
#endif

#ifndef rtmGetCtrlRateMdlRefTimingPtr
# define rtmGetCtrlRateMdlRefTimingPtr(rtm) ()
#endif

#ifndef rtmSetCtrlRateMdlRefTimingPtr
# define rtmSetCtrlRateMdlRefTimingPtr(rtm, val) ()
#endif

#ifndef rtmGetCtrlRateNumTicksToNextHit
# define rtmGetCtrlRateNumTicksToNextHit(rtm) ()
#endif

#ifndef rtmSetCtrlRateNumTicksToNextHit
# define rtmSetCtrlRateNumTicksToNextHit(rtm, val) ()
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        ()
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)   ()
#endif

#ifndef rtmGetDefaultParam
# define rtmGetDefaultParam(rtm)       ((rtm)->defaultParam)
#endif

#ifndef rtmSetDefaultParam
# define rtmSetDefaultParam(rtm, val)  ((rtm)->defaultParam = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetDirectFeedThrough
# define rtmGetDirectFeedThrough(rtm)  ((rtm)->Sizes.sysDirFeedThru)
#endif

#ifndef rtmSetDirectFeedThrough
# define rtmSetDirectFeedThrough(rtm, val) ((rtm)->Sizes.sysDirFeedThru = (val))
#endif

#ifndef rtmGetErrorStatusFlag
# define rtmGetErrorStatusFlag(rtm)    ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatusFlag
# define rtmSetErrorStatusFlag(rtm, val) ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetFinalTime
# define rtmSetFinalTime(rtm, val)     ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetFirstInitCondFlag
# define rtmGetFirstInitCondFlag(rtm)  ((rtm)->Timing.firstInitCondFlag)
#endif

#ifndef rtmSetFirstInitCondFlag
# define rtmSetFirstInitCondFlag(rtm, val) ((rtm)->Timing.firstInitCondFlag = (val))
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetMdlRefGlobalTID
# define rtmGetMdlRefGlobalTID(rtm)    ()
#endif

#ifndef rtmSetMdlRefGlobalTID
# define rtmSetMdlRefGlobalTID(rtm, val) ()
#endif

#ifndef rtmGetMdlRefTriggerTID
# define rtmGetMdlRefTriggerTID(rtm)   ()
#endif

#ifndef rtmSetMdlRefTriggerTID
# define rtmSetMdlRefTriggerTID(rtm, val) ()
#endif

#ifndef rtmGetModelMappingInfo
# define rtmGetModelMappingInfo(rtm)   ((rtm)->SpecialInfo.mappingInfo)
#endif

#ifndef rtmSetModelMappingInfo
# define rtmSetModelMappingInfo(rtm, val) ((rtm)->SpecialInfo.mappingInfo = (val))
#endif

#ifndef rtmGetModelName
# define rtmGetModelName(rtm)          ((rtm)->modelName)
#endif

#ifndef rtmSetModelName
# define rtmSetModelName(rtm, val)     ((rtm)->modelName = (val))
#endif

#ifndef rtmGetNonInlinedSFcns
# define rtmGetNonInlinedSFcns(rtm)    ()
#endif

#ifndef rtmSetNonInlinedSFcns
# define rtmSetNonInlinedSFcns(rtm, val) ()
#endif

#ifndef rtmGetNumBlockIO
# define rtmGetNumBlockIO(rtm)         ((rtm)->Sizes.numBlockIO)
#endif

#ifndef rtmSetNumBlockIO
# define rtmSetNumBlockIO(rtm, val)    ((rtm)->Sizes.numBlockIO = (val))
#endif

#ifndef rtmGetNumBlockParams
# define rtmGetNumBlockParams(rtm)     ((rtm)->Sizes.numBlockPrms)
#endif

#ifndef rtmSetNumBlockParams
# define rtmSetNumBlockParams(rtm, val) ((rtm)->Sizes.numBlockPrms = (val))
#endif

#ifndef rtmGetNumBlocks
# define rtmGetNumBlocks(rtm)          ((rtm)->Sizes.numBlocks)
#endif

#ifndef rtmSetNumBlocks
# define rtmSetNumBlocks(rtm, val)     ((rtm)->Sizes.numBlocks = (val))
#endif

#ifndef rtmGetNumContStates
# define rtmGetNumContStates(rtm)      ((rtm)->Sizes.numContStates)
#endif

#ifndef rtmSetNumContStates
# define rtmSetNumContStates(rtm, val) ((rtm)->Sizes.numContStates = (val))
#endif

#ifndef rtmGetNumDWork
# define rtmGetNumDWork(rtm)           ((rtm)->Sizes.numDwork)
#endif

#ifndef rtmSetNumDWork
# define rtmSetNumDWork(rtm, val)      ((rtm)->Sizes.numDwork = (val))
#endif

#ifndef rtmGetNumInputPorts
# define rtmGetNumInputPorts(rtm)      ((rtm)->Sizes.numIports)
#endif

#ifndef rtmSetNumInputPorts
# define rtmSetNumInputPorts(rtm, val) ((rtm)->Sizes.numIports = (val))
#endif

#ifndef rtmGetNumNonSampledZCs
# define rtmGetNumNonSampledZCs(rtm)   ((rtm)->Sizes.numNonSampZCs)
#endif

#ifndef rtmSetNumNonSampledZCs
# define rtmSetNumNonSampledZCs(rtm, val) ((rtm)->Sizes.numNonSampZCs = (val))
#endif

#ifndef rtmGetNumOutputPorts
# define rtmGetNumOutputPorts(rtm)     ((rtm)->Sizes.numOports)
#endif

#ifndef rtmSetNumOutputPorts
# define rtmSetNumOutputPorts(rtm, val) ((rtm)->Sizes.numOports = (val))
#endif

#ifndef rtmGetNumPeriodicContStates
# define rtmGetNumPeriodicContStates(rtm) ((rtm)->Sizes.numPeriodicContStates)
#endif

#ifndef rtmSetNumPeriodicContStates
# define rtmSetNumPeriodicContStates(rtm, val) ((rtm)->Sizes.numPeriodicContStates = (val))
#endif

#ifndef rtmGetNumSFcnParams
# define rtmGetNumSFcnParams(rtm)      ((rtm)->Sizes.numSFcnPrms)
#endif

#ifndef rtmSetNumSFcnParams
# define rtmSetNumSFcnParams(rtm, val) ((rtm)->Sizes.numSFcnPrms = (val))
#endif

#ifndef rtmGetNumSFunctions
# define rtmGetNumSFunctions(rtm)      ((rtm)->Sizes.numSFcns)
#endif

#ifndef rtmSetNumSFunctions
# define rtmSetNumSFunctions(rtm, val) ((rtm)->Sizes.numSFcns = (val))
#endif

#ifndef rtmGetNumSampleTimes
# define rtmGetNumSampleTimes(rtm)     ((rtm)->Sizes.numSampTimes)
#endif

#ifndef rtmSetNumSampleTimes
# define rtmSetNumSampleTimes(rtm, val) ((rtm)->Sizes.numSampTimes = (val))
#endif

#ifndef rtmGetNumU
# define rtmGetNumU(rtm)               ((rtm)->Sizes.numU)
#endif

#ifndef rtmSetNumU
# define rtmSetNumU(rtm, val)          ((rtm)->Sizes.numU = (val))
#endif

#ifndef rtmGetNumY
# define rtmGetNumY(rtm)               ((rtm)->Sizes.numY)
#endif

#ifndef rtmSetNumY
# define rtmSetNumY(rtm, val)          ((rtm)->Sizes.numY = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ()
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ()
#endif

#ifndef rtmGetOffsetTimeArray
# define rtmGetOffsetTimeArray(rtm)    ((rtm)->Timing.offsetTimesArray)
#endif

#ifndef rtmSetOffsetTimeArray
# define rtmSetOffsetTimeArray(rtm, val) ((rtm)->Timing.offsetTimesArray = (val))
#endif

#ifndef rtmGetOffsetTimePtr
# define rtmGetOffsetTimePtr(rtm)      ((rtm)->Timing.offsetTimes)
#endif

#ifndef rtmSetOffsetTimePtr
# define rtmSetOffsetTimePtr(rtm, val) ((rtm)->Timing.offsetTimes = (val))
#endif

#ifndef rtmGetOptions
# define rtmGetOptions(rtm)            ((rtm)->Sizes.options)
#endif

#ifndef rtmSetOptions
# define rtmSetOptions(rtm, val)       ((rtm)->Sizes.options = (val))
#endif

#ifndef rtmGetParamIsMalloced
# define rtmGetParamIsMalloced(rtm)    ()
#endif

#ifndef rtmSetParamIsMalloced
# define rtmSetParamIsMalloced(rtm, val) ()
#endif

#ifndef rtmGetPath
# define rtmGetPath(rtm)               ((rtm)->path)
#endif

#ifndef rtmSetPath
# define rtmSetPath(rtm, val)          ((rtm)->path = (val))
#endif

#ifndef rtmGetPerTaskSampleHits
# define rtmGetPerTaskSampleHits(rtm)  ((rtm)->Timing.RateInteraction)
#endif

#ifndef rtmSetPerTaskSampleHits
# define rtmSetPerTaskSampleHits(rtm, val) ((rtm)->Timing.RateInteraction = (val))
#endif

#ifndef rtmGetPerTaskSampleHitsArray
# define rtmGetPerTaskSampleHitsArray(rtm) ((rtm)->Timing.perTaskSampleHitsArray)
#endif

#ifndef rtmSetPerTaskSampleHitsArray
# define rtmSetPerTaskSampleHitsArray(rtm, val) ((rtm)->Timing.perTaskSampleHitsArray = (val))
#endif

#ifndef rtmGetPerTaskSampleHitsPtr
# define rtmGetPerTaskSampleHitsPtr(rtm) ((rtm)->Timing.perTaskSampleHits)
#endif

#ifndef rtmSetPerTaskSampleHitsPtr
# define rtmSetPerTaskSampleHitsPtr(rtm, val) ((rtm)->Timing.perTaskSampleHits = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetPrevZCSigState
# define rtmGetPrevZCSigState(rtm)     ((rtm)->prevZCSigState)
#endif

#ifndef rtmSetPrevZCSigState
# define rtmSetPrevZCSigState(rtm, val) ((rtm)->prevZCSigState = (val))
#endif

#ifndef rtmGetRTWExtModeInfo
# define rtmGetRTWExtModeInfo(rtm)     ((rtm)->extModeInfo)
#endif

#ifndef rtmSetRTWExtModeInfo
# define rtmSetRTWExtModeInfo(rtm, val) ((rtm)->extModeInfo = (val))
#endif

#ifndef rtmGetRTWGeneratedSFcn
# define rtmGetRTWGeneratedSFcn(rtm)   ((rtm)->Sizes.rtwGenSfcn)
#endif

#ifndef rtmSetRTWGeneratedSFcn
# define rtmSetRTWGeneratedSFcn(rtm, val) ((rtm)->Sizes.rtwGenSfcn = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmSetRTWLogInfo
# define rtmSetRTWLogInfo(rtm, val)    ((rtm)->rtwLogInfo = (val))
#endif

#ifndef rtmGetRTWRTModelMethodsInfo
# define rtmGetRTWRTModelMethodsInfo(rtm) ()
#endif

#ifndef rtmSetRTWRTModelMethodsInfo
# define rtmSetRTWRTModelMethodsInfo(rtm, val) ()
#endif

#ifndef rtmGetRTWSfcnInfo
# define rtmGetRTWSfcnInfo(rtm)        ((rtm)->sfcnInfo)
#endif

#ifndef rtmSetRTWSfcnInfo
# define rtmSetRTWSfcnInfo(rtm, val)   ((rtm)->sfcnInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfo
# define rtmGetRTWSolverInfo(rtm)      ((rtm)->solverInfo)
#endif

#ifndef rtmSetRTWSolverInfo
# define rtmSetRTWSolverInfo(rtm, val) ((rtm)->solverInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfoPtr
# define rtmGetRTWSolverInfoPtr(rtm)   ((rtm)->solverInfoPtr)
#endif

#ifndef rtmSetRTWSolverInfoPtr
# define rtmSetRTWSolverInfoPtr(rtm, val) ((rtm)->solverInfoPtr = (val))
#endif

#ifndef rtmGetReservedForXPC
# define rtmGetReservedForXPC(rtm)     ((rtm)->SpecialInfo.xpcData)
#endif

#ifndef rtmSetReservedForXPC
# define rtmSetReservedForXPC(rtm, val) ((rtm)->SpecialInfo.xpcData = (val))
#endif

#ifndef rtmGetRootDWork
# define rtmGetRootDWork(rtm)          ((rtm)->dwork)
#endif

#ifndef rtmSetRootDWork
# define rtmSetRootDWork(rtm, val)     ((rtm)->dwork = (val))
#endif

#ifndef rtmGetSFunctions
# define rtmGetSFunctions(rtm)         ((rtm)->childSfunctions)
#endif

#ifndef rtmSetSFunctions
# define rtmSetSFunctions(rtm, val)    ((rtm)->childSfunctions = (val))
#endif

#ifndef rtmGetSampleHitArray
# define rtmGetSampleHitArray(rtm)     ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmSetSampleHitArray
# define rtmSetSampleHitArray(rtm, val) ((rtm)->Timing.sampleHitArray = (val))
#endif

#ifndef rtmGetSampleHitPtr
# define rtmGetSampleHitPtr(rtm)       ((rtm)->Timing.sampleHits)
#endif

#ifndef rtmSetSampleHitPtr
# define rtmSetSampleHitPtr(rtm, val)  ((rtm)->Timing.sampleHits = (val))
#endif

#ifndef rtmGetSampleTimeArray
# define rtmGetSampleTimeArray(rtm)    ((rtm)->Timing.sampleTimesArray)
#endif

#ifndef rtmSetSampleTimeArray
# define rtmSetSampleTimeArray(rtm, val) ((rtm)->Timing.sampleTimesArray = (val))
#endif

#ifndef rtmGetSampleTimePtr
# define rtmGetSampleTimePtr(rtm)      ((rtm)->Timing.sampleTimes)
#endif

#ifndef rtmSetSampleTimePtr
# define rtmSetSampleTimePtr(rtm, val) ((rtm)->Timing.sampleTimes = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDArray
# define rtmGetSampleTimeTaskIDArray(rtm) ((rtm)->Timing.sampleTimeTaskIDArray)
#endif

#ifndef rtmSetSampleTimeTaskIDArray
# define rtmSetSampleTimeTaskIDArray(rtm, val) ((rtm)->Timing.sampleTimeTaskIDArray = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDPtr
# define rtmGetSampleTimeTaskIDPtr(rtm) ((rtm)->Timing.sampleTimeTaskIDPtr)
#endif

#ifndef rtmSetSampleTimeTaskIDPtr
# define rtmSetSampleTimeTaskIDPtr(rtm, val) ((rtm)->Timing.sampleTimeTaskIDPtr = (val))
#endif

#ifndef rtmGetSelf
# define rtmGetSelf(rtm)               ()
#endif

#ifndef rtmSetSelf
# define rtmSetSelf(rtm, val)          ()
#endif

#ifndef rtmGetSimMode
# define rtmGetSimMode(rtm)            ((rtm)->simMode)
#endif

#ifndef rtmSetSimMode
# define rtmSetSimMode(rtm, val)       ((rtm)->simMode = (val))
#endif

#ifndef rtmGetSimTimeStep
# define rtmGetSimTimeStep(rtm)        ((rtm)->Timing.simTimeStep)
#endif

#ifndef rtmSetSimTimeStep
# define rtmSetSimTimeStep(rtm, val)   ((rtm)->Timing.simTimeStep = (val))
#endif

#ifndef rtmGetStartTime
# define rtmGetStartTime(rtm)          ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetStartTime
# define rtmSetStartTime(rtm, val)     ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetStepSize
# define rtmGetStepSize(rtm)           ((rtm)->Timing.stepSize)
#endif

#ifndef rtmSetStepSize
# define rtmSetStepSize(rtm, val)      ((rtm)->Timing.stepSize = (val))
#endif

#ifndef rtmGetStopRequestedFlag
# define rtmGetStopRequestedFlag(rtm)  ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequestedFlag
# define rtmSetStopRequestedFlag(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetTaskCounters
# define rtmGetTaskCounters(rtm)       ((rtm)->Timing.TaskCounters)
#endif

#ifndef rtmSetTaskCounters
# define rtmSetTaskCounters(rtm, val)  ((rtm)->Timing.TaskCounters = (val))
#endif

#ifndef rtmGetTaskTimeArray
# define rtmGetTaskTimeArray(rtm)      ((rtm)->Timing.tArray)
#endif

#ifndef rtmSetTaskTimeArray
# define rtmSetTaskTimeArray(rtm, val) ((rtm)->Timing.tArray = (val))
#endif

#ifndef rtmGetTimePtr
# define rtmGetTimePtr(rtm)            ((rtm)->Timing.t)
#endif

#ifndef rtmSetTimePtr
# define rtmSetTimePtr(rtm, val)       ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTimingData
# define rtmGetTimingData(rtm)         ((rtm)->Timing.timingData)
#endif

#ifndef rtmSetTimingData
# define rtmSetTimingData(rtm, val)    ((rtm)->Timing.timingData = (val))
#endif

#ifndef rtmGetU
# define rtmGetU(rtm)                  ((rtm)->inputs)
#endif

#ifndef rtmSetU
# define rtmSetU(rtm, val)             ((rtm)->inputs = (val))
#endif

#ifndef rtmGetVarNextHitTimesListPtr
# define rtmGetVarNextHitTimesListPtr(rtm) ((rtm)->Timing.varNextHitTimesList)
#endif

#ifndef rtmSetVarNextHitTimesListPtr
# define rtmSetVarNextHitTimesListPtr(rtm, val) ((rtm)->Timing.varNextHitTimesList = (val))
#endif

#ifndef rtmGetY
# define rtmGetY(rtm)                  ((rtm)->outputs)
#endif

#ifndef rtmSetY
# define rtmSetY(rtm, val)             ((rtm)->outputs = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetZCSignalValues
# define rtmGetZCSignalValues(rtm)     ((rtm)->zcSignalValues)
#endif

#ifndef rtmSetZCSignalValues
# define rtmSetZCSignalValues(rtm, val) ((rtm)->zcSignalValues = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
# define rtmGet_TimeOfLastOutput(rtm)  ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmSet_TimeOfLastOutput
# define rtmSet_TimeOfLastOutput(rtm, val) ((rtm)->Timing.timeOfLastOutput = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGettimingBridge
# define rtmGettimingBridge(rtm)       ()
#endif

#ifndef rtmSettimingBridge
# define rtmSettimingBridge(rtm, val)  ()
#endif

#ifndef rtmGetChecksumVal
# define rtmGetChecksumVal(rtm, idx)   ((rtm)->Sizes.checksums[idx])
#endif

#ifndef rtmSetChecksumVal
# define rtmSetChecksumVal(rtm, idx, val) ((rtm)->Sizes.checksums[idx] = (val))
#endif

#ifndef rtmGetDWork
# define rtmGetDWork(rtm, idx)         ((rtm)->dwork[idx])
#endif

#ifndef rtmSetDWork
# define rtmSetDWork(rtm, idx, val)    ((rtm)->dwork[idx] = (val))
#endif

#ifndef rtmGetOffsetTime
# define rtmGetOffsetTime(rtm, idx)    ((rtm)->Timing.offsetTimes[idx])
#endif

#ifndef rtmSetOffsetTime
# define rtmSetOffsetTime(rtm, idx, val) ((rtm)->Timing.offsetTimes[idx] = (val))
#endif

#ifndef rtmGetSFunction
# define rtmGetSFunction(rtm, idx)     ((rtm)->childSfunctions[idx])
#endif

#ifndef rtmSetSFunction
# define rtmSetSFunction(rtm, idx, val) ((rtm)->childSfunctions[idx] = (val))
#endif

#ifndef rtmGetSampleTime
# define rtmGetSampleTime(rtm, idx)    ((rtm)->Timing.sampleTimes[idx])
#endif

#ifndef rtmSetSampleTime
# define rtmSetSampleTime(rtm, idx, val) ((rtm)->Timing.sampleTimes[idx] = (val))
#endif

#ifndef rtmGetSampleTimeTaskID
# define rtmGetSampleTimeTaskID(rtm, idx) ((rtm)->Timing.sampleTimeTaskIDPtr[idx])
#endif

#ifndef rtmSetSampleTimeTaskID
# define rtmSetSampleTimeTaskID(rtm, idx, val) ((rtm)->Timing.sampleTimeTaskIDPtr[idx] = (val))
#endif

#ifndef rtmGetVarNextHitTimeList
# define rtmGetVarNextHitTimeList(rtm, idx) ((rtm)->Timing.varNextHitTimesList[idx])
#endif

#ifndef rtmSetVarNextHitTimeList
# define rtmSetVarNextHitTimeList(rtm, idx, val) ((rtm)->Timing.varNextHitTimesList[idx] = (val))
#endif

#ifndef rtmIsContinuousTask
# define rtmIsContinuousTask(rtm, tid) ((tid) <= 1)
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmSetFirstInitCond
# define rtmSetFirstInitCond(rtm, val) ((rtm)->Timing.firstInitCondFlag = (val))
#endif

#ifndef rtmIsFirstInitCond
# define rtmIsFirstInitCond(rtm)       ((rtm)->Timing.firstInitCondFlag)
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmIsSampleHit
# define rtmIsSampleHit(rtm, sti, tid) (((rtm)->Timing.sampleTimeTaskIDPtr[sti] == (tid)))
#endif

#ifndef rtmStepTask
# define rtmStepTask(rtm, idx)         ((rtm)->Timing.TaskCounters.TID[(idx)] == 0)
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmSetT
# define rtmSetT(rtm, val)                                       /* Do Nothing */
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetTFinal
# define rtmSetTFinal(rtm, val)        ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTStart
# define rtmGetTStart(rtm)             ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetTStart
# define rtmSetTStart(rtm, val)        ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmTaskCounter
# define rtmTaskCounter(rtm, idx)      ((rtm)->Timing.TaskCounters.TID[(idx)])
#endif

#ifndef rtmGetTaskTime
# define rtmGetTaskTime(rtm, sti)      (rtmGetTPtr((rtm))[(rtm)->Timing.sampleTimeTaskIDPtr[sti]])
#endif

#ifndef rtmSetTaskTime
# define rtmSetTaskTime(rtm, sti, val) (rtmGetTPtr((rtm))[sti] = (val))
#endif

#ifndef rtmGetTimeOfLastOutput
# define rtmGetTimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifdef rtmGetRTWSolverInfo
#undef rtmGetRTWSolverInfo
#endif

#define rtmGetRTWSolverInfo(rtm)       &((rtm)->solverInfo)

/* Definition for use in the target main file */
#define ARDroneTTSim_rtModel           RT_MODEL_ARDroneTTSim_T

/* Block signals for system '<S1>/normalize angle  between -pi and pi radians' */
typedef struct {
  real_T angleOut;      /* '<S1>/normalize angle  between -pi and pi radians' */
} B_normalizeanglebetweenpiandp_T;

/* Block signals (default storage) */
typedef struct {
  real_T gain1;                        /* '<S1>/gain1' */
  real_T deg2rad1;                     /* '<S1>/deg 2 rad1' */
  real_T gain;                         /* '<S1>/gain' */
  real_T deg2rad;                      /* '<S1>/deg 2 rad' */
  real_T deg2rad2;                     /* '<S1>/deg 2 rad2' */
  real_T deg2rad_p[3];                 /* '<Root>/deg 2 rad' */
  real_T StateSpace3;                  /* '<S1>/State-Space3' */
  real_T StateSpace2;                  /* '<S1>/State-Space2' */
  real_T Constant1[2];                 /* '<S19>/Constant1' */
  real_T Integrator[2];                /* '<S19>/Integrator' */
  real_T StateSpace5;                  /* '<S1>/State-Space5' */
  real_T Time;                         /* '<Root>/Time' */
  real_T deg2rad1_l;                   /* '<S5>/deg 2 rad1' */
  real_T RTrans4;                      /* '<Root>/RTrans4' */
  real_T Saturation1;                  /* '<S1>/Saturation 1' */
  real_T RTrans3;                      /* '<Root>/RTrans3' */
  real_T Saturation2;                  /* '<S1>/Saturation 2' */
  real_T RTrans2;                      /* '<Root>/RTrans2' */
  real_T Saturation3;                  /* '<S1>/Saturation 3' */
  real_T RTrans1;                      /* '<Root>/RTrans1' */
  real_T Saturation4;                  /* '<S1>/Saturation 4' */
  real_T totalcommunicationtimedelay[4];
                                     /* '<S1>/total communication time delay' */
  real_T ManualSwitch3[4];             /* '<S3>/Manual Switch3' */
  real_T Vel_xy[2];
          /* '<S19>/Velocity from vehicle body frame  to inertial NED  frame' */
  real_T dot_xi[4];                    /* '<S3>/MATLAB Function' */
  boolean_T Compare;                   /* '<S8>/Compare' */
  B_normalizeanglebetweenpiandp_T sf_normalizeanglebetweenpiand_c;
                        /* '<S9>/normalize angle  between -pi and pi radians' */
  B_normalizeanglebetweenpiandp_T sf_normalizeanglebetweenpiand_e;
                       /* '<S12>/normalize angle  between -pi and pi radians' */
  B_normalizeanglebetweenpiandp_T sf_normalizeanglebetweenpiandpi;
                        /* '<S1>/normalize angle  between -pi and pi radians' */
} B_ARDroneTTSim_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T DiscreteTimeIntegrator_DSTATE[4];/* '<S3>/Discrete-Time Integrator' */
  real_T RTrans5_1_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans5_2_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans5_3_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans5_4_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans5_5_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans5_6_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans5_7_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans5_8_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans5_9_Buffer;             /* '<Root>/RTrans5' */
  real_T RTrans7_Buffer;               /* '<Root>/RTrans7' */
  real_T RTrans4_Buffer0;              /* '<Root>/RTrans4' */
  real_T RTrans3_Buffer0;              /* '<Root>/RTrans3' */
  real_T RTrans2_Buffer0;              /* '<Root>/RTrans2' */
  real_T RTrans1_Buffer0;              /* '<Root>/RTrans1' */
  real_T t0;                           /* '<Root>/MATLAB Function1' */
  real_T p0[3];                        /* '<Root>/MATLAB Function1' */
  real_T previous_status;              /* '<Root>/MATLAB Function1' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[32768];
  } totalcommunicationtimedelay_RWO; /* '<S1>/total communication time delay' */

  struct {
    void *LoggedData;
  } Euleranglesdeg_PWORK;              /* '<Root>/Euler angles (deg)' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK;                 /* '<Root>/To Workspace' */

  struct {
    void *LoggedData;
  } Heightm_PWORK;                     /* '<S5>/Height (m)' */

  struct {
    void *LoggedData;
  } InertialpotitionalongXem_PWORK;  /* '<S5>/Inertial potition along Xe (m)' */

  struct {
    void *LoggedData;
  } InertialpotitionalongYem_PWORK;  /* '<S5>/Inertial potition along Ye (m)' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK_d;               /* '<S5>/To Workspace' */

  struct {
    void *LoggedData;
  } headingdeg_PWORK;                  /* '<S5>/heading (deg)' */

  struct {
    void *LoggedData;
  } ToWorkspace1_PWORK;                /* '<S5>/To Workspace1' */

  struct {
    void *LoggedData[2];
  } pdm_PWORK;                         /* '<Root>/pd (m)' */

  struct {
    void *TUbufferPtrs[8];
  } totalcommunicationtimedelay_PWO; /* '<S1>/total communication time delay' */

  struct {
    void *LoggedData;
  } error_phi_PWORK;                   /* '<S9>/error_phi' */

  struct {
    void *LoggedData;
  } error_psi_PWORK;                   /* '<S9>/error_psi' */

  struct {
    void *LoggedData;
  } error_theta_PWORK;                 /* '<S9>/error_theta' */

  struct {
    void *LoggedData;
  } errorx_PWORK;                      /* '<S9>/errorx' */

  struct {
    void *LoggedData;
  } errory_PWORK;                      /* '<S9>/errory' */

  struct {
    void *LoggedData;
  } errorz_PWORK;                      /* '<S9>/errorz' */

  struct {
    void *LoggedData;
  } psipsi_d_PWORK;                    /* '<S9>/psi  psi_d' */

  int_T Integrator_IWORK;              /* '<S19>/Integrator' */
  struct {
    int_T Tail[4];
    int_T Head[4];
    int_T Last[4];
    int_T CircularBufSize[4];
  } totalcommunicationtimedelay_IWO; /* '<S1>/total communication time delay' */

  int8_T DiscreteTimeIntegrator_PrevRese;/* '<S3>/Discrete-Time Integrator' */
  boolean_T t0_not_empty;              /* '<Root>/MATLAB Function1' */
  boolean_T previous_status_not_empty; /* '<Root>/MATLAB Function1' */
} DW_ARDroneTTSim_T;

/* Continuous states (default storage) */
typedef struct {
  real_T StateSpace_CSTATE[2];         /* '<S1>/State-Space' */
  real_T StateSpace1_CSTATE[2];        /* '<S1>/State-Space1' */
  real_T StateSpace4_CSTATE;           /* '<S1>/State-Space4' */
  real_T StateSpace3_CSTATE;           /* '<S1>/State-Space3' */
  real_T StateSpace2_CSTATE;           /* '<S1>/State-Space2' */
  real_T Integrator_CSTATE[2];         /* '<S19>/Integrator' */
  real_T StateSpace5_CSTATE[2];        /* '<S1>/State-Space5' */
} X_ARDroneTTSim_T;

/* State derivatives (default storage) */
typedef struct {
  real_T StateSpace_CSTATE[2];         /* '<S1>/State-Space' */
  real_T StateSpace1_CSTATE[2];        /* '<S1>/State-Space1' */
  real_T StateSpace4_CSTATE;           /* '<S1>/State-Space4' */
  real_T StateSpace3_CSTATE;           /* '<S1>/State-Space3' */
  real_T StateSpace2_CSTATE;           /* '<S1>/State-Space2' */
  real_T Integrator_CSTATE[2];         /* '<S19>/Integrator' */
  real_T StateSpace5_CSTATE[2];        /* '<S1>/State-Space5' */
} XDot_ARDroneTTSim_T;

/* State disabled  */
typedef struct {
  boolean_T StateSpace_CSTATE[2];      /* '<S1>/State-Space' */
  boolean_T StateSpace1_CSTATE[2];     /* '<S1>/State-Space1' */
  boolean_T StateSpace4_CSTATE;        /* '<S1>/State-Space4' */
  boolean_T StateSpace3_CSTATE;        /* '<S1>/State-Space3' */
  boolean_T StateSpace2_CSTATE;        /* '<S1>/State-Space2' */
  boolean_T Integrator_CSTATE[2];      /* '<S19>/Integrator' */
  boolean_T StateSpace5_CSTATE[2];     /* '<S1>/State-Space5' */
} XDis_ARDroneTTSim_T;

#ifndef ODE1_INTG
#define ODE1_INTG

/* ODE1 Integration Data */
typedef struct {
  real_T *f[1];                        /* derivatives */
} ODE1_IntgData;

#endif

/* Backward compatible GRT Identifiers */
#define rtB                            ARDroneTTSim_B
#define BlockIO                        B_ARDroneTTSim_T
#define rtX                            ARDroneTTSim_X
#define ContinuousStates               X_ARDroneTTSim_T
#define rtXdot                         ARDroneTTSim_XDot
#define StateDerivatives               XDot_ARDroneTTSim_T
#define tXdis                          ARDroneTTSim_XDis
#define StateDisabled                  XDis_ARDroneTTSim_T
#define rtP                            ARDroneTTSim_P
#define Parameters                     P_ARDroneTTSim_T
#define rtDWork                        ARDroneTTSim_DW
#define D_Work                         DW_ARDroneTTSim_T

/* Parameters (default storage) */
struct P_ARDroneTTSim_T_ {
  real_T K[18];                        /* Variable: K
                                        * Referenced by: '<S3>/Constant1'
                                        */
  real_T k_w;                          /* Variable: k_w
                                        * Referenced by: '<S3>/Constant'
                                        */
  real_T timeDelay;                    /* Variable: timeDelay
                                        * Referenced by: '<S1>/time delay'
                                        */
  real_T commandsenabled_Value;        /* Expression: 1
                                        * Referenced by: '<Root>/commands enabled'
                                        */
  real_T commandsdisabled_Value;       /* Expression: 0
                                        * Referenced by: '<Root>/commands disabled'
                                        */
  real_T Gain2_Gain;                   /* Expression: 1
                                        * Referenced by: '<S15>/Gain2'
                                        */
  real_T Gain_Gain;                    /* Expression: -0.5
                                        * Referenced by: '<S11>/Gain'
                                        */
  real_T proportionalcontrolgainyaw_Gain;/* Expression: 1.5
                                          * Referenced by: '<S12>/proportional control gain - yaw'
                                          */
  real_T proportionalcontrolgain_Gain; /* Expression: 1
                                        * Referenced by: '<S13>/proportional control gain'
                                        */
  real_T Gain3_Gain;                   /* Expression: 1
                                        * Referenced by: '<S15>/Gain3'
                                        */
  real_T Gain1_Gain;                   /* Expression: 0.4
                                        * Referenced by: '<S14>/Gain1'
                                        */
  real_T StateSpace_A[3];              /* Computed Parameter: StateSpace_A
                                        * Referenced by: '<S1>/State-Space'
                                        */
  real_T StateSpace_B;                 /* Computed Parameter: StateSpace_B
                                        * Referenced by: '<S1>/State-Space'
                                        */
  real_T StateSpace_C[2];              /* Computed Parameter: StateSpace_C
                                        * Referenced by: '<S1>/State-Space'
                                        */
  real_T StateSpace_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S1>/State-Space'
                                        */
  real_T gain1_Gain;                   /* Expression: 1
                                        * Referenced by: '<S1>/gain1'
                                        */
  real_T deg2rad1_Gain;                /* Expression: 180/pi
                                        * Referenced by: '<S1>/deg 2 rad1'
                                        */
  real_T StateSpace1_A[3];             /* Computed Parameter: StateSpace1_A
                                        * Referenced by: '<S1>/State-Space1'
                                        */
  real_T StateSpace1_B;                /* Computed Parameter: StateSpace1_B
                                        * Referenced by: '<S1>/State-Space1'
                                        */
  real_T StateSpace1_C[2];             /* Computed Parameter: StateSpace1_C
                                        * Referenced by: '<S1>/State-Space1'
                                        */
  real_T StateSpace1_InitialCondition; /* Expression: 0
                                        * Referenced by: '<S1>/State-Space1'
                                        */
  real_T gain_Gain;                    /* Expression: 1
                                        * Referenced by: '<S1>/gain'
                                        */
  real_T deg2rad_Gain;                 /* Expression: 180/pi
                                        * Referenced by: '<S1>/deg 2 rad'
                                        */
  real_T StateSpace4_A;                /* Computed Parameter: StateSpace4_A
                                        * Referenced by: '<S1>/State-Space4'
                                        */
  real_T StateSpace4_B;                /* Computed Parameter: StateSpace4_B
                                        * Referenced by: '<S1>/State-Space4'
                                        */
  real_T StateSpace4_C;                /* Computed Parameter: StateSpace4_C
                                        * Referenced by: '<S1>/State-Space4'
                                        */
  real_T StateSpace4_InitialCondition; /* Expression: 0
                                        * Referenced by: '<S1>/State-Space4'
                                        */
  real_T deg2rad2_Gain;                /* Expression: 180/pi
                                        * Referenced by: '<S1>/deg 2 rad2'
                                        */
  real_T deg2rad_Gain_a;               /* Expression: pi/180
                                        * Referenced by: '<Root>/deg 2 rad'
                                        */
  real_T StateSpace3_A;                /* Computed Parameter: StateSpace3_A
                                        * Referenced by: '<S1>/State-Space3'
                                        */
  real_T StateSpace3_B;                /* Computed Parameter: StateSpace3_B
                                        * Referenced by: '<S1>/State-Space3'
                                        */
  real_T StateSpace3_C;                /* Computed Parameter: StateSpace3_C
                                        * Referenced by: '<S1>/State-Space3'
                                        */
  real_T StateSpace3_InitialCondition; /* Expression: 0
                                        * Referenced by: '<S1>/State-Space3'
                                        */
  real_T StateSpace2_A;                /* Computed Parameter: StateSpace2_A
                                        * Referenced by: '<S1>/State-Space2'
                                        */
  real_T StateSpace2_B;                /* Computed Parameter: StateSpace2_B
                                        * Referenced by: '<S1>/State-Space2'
                                        */
  real_T StateSpace2_C;                /* Computed Parameter: StateSpace2_C
                                        * Referenced by: '<S1>/State-Space2'
                                        */
  real_T StateSpace2_InitialCondition; /* Expression: 0
                                        * Referenced by: '<S1>/State-Space2'
                                        */
  real_T TheARDronesendszeroforthevertic;/* Expression: 0
                                          * Referenced by: '<S1>/The ARDrone sends zero for the vertical velocity.  '
                                          */
  real_T Constant1_Value[2];           /* Expression: [0 0]
                                        * Referenced by: '<S19>/Constant1'
                                        */
  real_T StateSpace5_A[3];             /* Computed Parameter: StateSpace5_A
                                        * Referenced by: '<S1>/State-Space5'
                                        */
  real_T StateSpace5_B;                /* Computed Parameter: StateSpace5_B
                                        * Referenced by: '<S1>/State-Space5'
                                        */
  real_T StateSpace5_C[2];             /* Computed Parameter: StateSpace5_C
                                        * Referenced by: '<S1>/State-Space5'
                                        */
  real_T StateSpace5_InitialCondition; /* Expression: 0
                                        * Referenced by: '<S1>/State-Space5'
                                        */
  real_T deg2rad1_Gain_n;              /* Expression: 180/pi
                                        * Referenced by: '<S5>/deg 2 rad1'
                                        */
  real_T RTrans4_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<Root>/RTrans4'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 1
                                        * Referenced by: '<S1>/Saturation 1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: -1
                                        * Referenced by: '<S1>/Saturation 1'
                                        */
  real_T RTrans3_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<Root>/RTrans3'
                                        */
  real_T Saturation2_UpperSat;         /* Expression: 1
                                        * Referenced by: '<S1>/Saturation 2'
                                        */
  real_T Saturation2_LowerSat;         /* Expression: -1
                                        * Referenced by: '<S1>/Saturation 2'
                                        */
  real_T RTrans2_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<Root>/RTrans2'
                                        */
  real_T Saturation3_UpperSat;         /* Expression: 1
                                        * Referenced by: '<S1>/Saturation 3'
                                        */
  real_T Saturation3_LowerSat;         /* Expression: -1
                                        * Referenced by: '<S1>/Saturation 3'
                                        */
  real_T RTrans1_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<Root>/RTrans1'
                                        */
  real_T Saturation4_UpperSat;         /* Expression: 1
                                        * Referenced by: '<S1>/Saturation 4'
                                        */
  real_T Saturation4_LowerSat;         /* Expression: -1
                                        * Referenced by: '<S1>/Saturation 4'
                                        */
  real_T totalcommunicationtimedelay_Max;/* Expression: 10
                                          * Referenced by: '<S1>/total communication time delay'
                                          */
  real_T totalcommunicationtimedelay_Ini[4];/* Expression: [0 0 0 0]
                                             * Referenced by: '<S1>/total communication time delay'
                                             */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<S8>/Constant'
                                        */
  real_T DiscreteTimeIntegrator_gainval;
                           /* Computed Parameter: DiscreteTimeIntegrator_gainval
                            * Referenced by: '<S3>/Discrete-Time Integrator'
                            */
  real_T DiscreteTimeIntegrator_IC[4]; /* Expression: [0;0;0;0]
                                        * Referenced by: '<S3>/Discrete-Time Integrator'
                                        */
  uint8_T ManualSwitch2_CurrentSetting;
                             /* Computed Parameter: ManualSwitch2_CurrentSetting
                              * Referenced by: '<Root>/Manual Switch2'
                              */
  uint8_T ManualSwitch3_CurrentSetting;
                             /* Computed Parameter: ManualSwitch3_CurrentSetting
                              * Referenced by: '<S3>/Manual Switch3'
                              */
};

/* Real-time Model Data Structure */
struct tag_RTM_ARDroneTTSim_T {
  const char_T *path;
  const char_T *modelName;
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWLogInfo *rtwLogInfo;
  RTWExtModeInfo *extModeInfo;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;
  void *blockIO;
  const void *constBlockIO;
  void *defaultParam;
  ZCSigState *prevZCSigState;
  real_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  void *zcSignalValues;
  void *inputs;
  void *outputs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeF[1][11];
  ODE1_IntgData intgData;
  void *dwork;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T checksums[4];
    uint32_T options;
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * SpecialInfo:
   * The following substructure contains special information
   * related to other components that are dependent on RTW.
   */
  struct {
    const void *mappingInfo;
    void *xpcData;
  } SpecialInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T stepSize1;
    uint32_T clockTick2;
    uint32_T clockTickH2;
    time_T stepSize2;
    boolean_T firstInitCondFlag;
    struct {
      uint8_T TID[3];
    } TaskCounters;

    struct {
      boolean_T TID1_2;
    } RateInteraction;

    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    void *timingData;
    real_T *varNextHitTimesList;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[3];
    time_T offsetTimesArray[3];
    int_T sampleTimeTaskIDArray[3];
    int_T sampleHitArray[3];
    int_T perTaskSampleHitsArray[9];
    time_T tArray[3];
  } Timing;
};

/* Block parameters (default storage) */
extern P_ARDroneTTSim_T ARDroneTTSim_P;

/* Block signals (default storage) */
extern B_ARDroneTTSim_T ARDroneTTSim_B;

/* Continuous states (default storage) */
extern X_ARDroneTTSim_T ARDroneTTSim_X;

/* Block states (default storage) */
extern DW_ARDroneTTSim_T ARDroneTTSim_DW;

/* External function called from main */
extern time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
  ;

/* Model entry point functions */
extern void ARDroneTTSim_initialize(void);
extern void ARDroneTTSim_output0(void);
extern void ARDroneTTSim_update0(void);
extern void ARDroneTTSim_output(int_T tid);
extern void ARDroneTTSim_update(int_T tid);
extern void ARDroneTTSim_terminate(void);

/*====================*
 * External functions *
 *====================*/
extern ARDroneTTSim_rtModel *ARDroneTTSim(void);
extern void MdlInitializeSizes(void);
extern void MdlInitializeSampleTimes(void);
extern void MdlInitialize(void);
extern void MdlStart(void);
extern void MdlOutputs(int_T tid);
extern void MdlUpdate(int_T tid);
extern void MdlTerminate(void);

/* Real-time Model object */
extern RT_MODEL_ARDroneTTSim_T *const ARDroneTTSim_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'ARDroneTTSim'
 * '<S1>'   : 'ARDroneTTSim/ARDrone Simulation Block'
 * '<S2>'   : 'ARDroneTTSim/MATLAB Function1'
 * '<S3>'   : 'ARDroneTTSim/Outer-loop Controller'
 * '<S4>'   : 'ARDroneTTSim/Position estimation Important:This block provides an  inaccurate extimation of position  based on  velocity information. '
 * '<S5>'   : 'ARDroneTTSim/Visualization of Drone states'
 * '<S6>'   : 'ARDroneTTSim/ARDrone Simulation Block/normalize angle  between -pi and pi radians'
 * '<S7>'   : 'ARDroneTTSim/Outer-loop Controller/Baseline Controller'
 * '<S8>'   : 'ARDroneTTSim/Outer-loop Controller/Compare To Zero'
 * '<S9>'   : 'ARDroneTTSim/Outer-loop Controller/Error Scopes'
 * '<S10>'  : 'ARDroneTTSim/Outer-loop Controller/MATLAB Function'
 * '<S11>'  : 'ARDroneTTSim/Outer-loop Controller/Baseline Controller/Forward velocity controller '
 * '<S12>'  : 'ARDroneTTSim/Outer-loop Controller/Baseline Controller/Heading controller'
 * '<S13>'  : 'ARDroneTTSim/Outer-loop Controller/Baseline Controller/Height controller '
 * '<S14>'  : 'ARDroneTTSim/Outer-loop Controller/Baseline Controller/Lateral velocity controller'
 * '<S15>'  : 'ARDroneTTSim/Outer-loop Controller/Baseline Controller/Position controller'
 * '<S16>'  : 'ARDroneTTSim/Outer-loop Controller/Baseline Controller/Heading controller/normalize angle  between -pi and pi radians'
 * '<S17>'  : 'ARDroneTTSim/Outer-loop Controller/Baseline Controller/Position controller/Coordinate trnasformation  from inertial frame to body frame '
 * '<S18>'  : 'ARDroneTTSim/Outer-loop Controller/Error Scopes/normalize angle  between -pi and pi radians'
 * '<S19>'  : 'ARDroneTTSim/Position estimation Important:This block provides an  inaccurate extimation of position  based on  velocity information. /Position estimation'
 * '<S20>'  : 'ARDroneTTSim/Position estimation Important:This block provides an  inaccurate extimation of position  based on  velocity information. /Position estimation/Velocity from vehicle body frame  to inertial NED  frame'
 */
#endif                                 /* RTW_HEADER_ARDroneTTSim_h_ */
