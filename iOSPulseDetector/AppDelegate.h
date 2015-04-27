/*****************************************************************************
 *   AppDelegate.h
 ******************************************************************************
 *   by Melissa Mozifian, 25th Dec 2014
 ******************************************************************************
 *
 *
 *  
 *
 *****************************************************************************/

#import <UIKit/UIKit.h>

@interface AppDelegate : UIResponder <UIApplicationDelegate>

@property (strong, nonatomic) UIWindow *window;

@property (nonatomic, strong) UIAlertController* alert;
@property (nonatomic, strong) UIAlertAction* allowCamAccess;
@property (nonatomic, strong) UIAlertAction* cancelAccess;

@end
