//
//  VFMap.h
//  VFMap
//
//  Copyright 2006-2021 All rights reserved.
//

#import "VF.h"
#import "VFStateSystem.h"
@class VFMapStimTable;

@interface VFMap:LLTaskPlugIn {

    NSPoint                 currentEyesUnits[kEyes];
    NSWindowController      *summaryController;
    NSArray                 *topLevelObjects;
}

@property (readonly, strong) LLDigitalOut *digitalOut;
@property (strong) NSArray *observerKeyArray;
@property (strong) VFStimuli *stimuli;
@property (readonly, strong) VFMapStimTable *mapStimTable;

- (IBAction)doJuice:(id)sender;
- (IBAction)doReset:(id)sender;

@end

