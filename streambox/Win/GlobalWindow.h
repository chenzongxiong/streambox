#ifndef GLOBALWINDOW_H
#define GLOBALWINDOW_H
/**
 * The default window into which all data is placed (via {@link GlobalWindows}).
 */
class GlobalWindow : BoundedWindow{
public:

	/**
	 * Singleton instance of {@link GlobalWindow}.
	 */
	//GlobalWindow(){}
        static GlobalWindow INSTANCE;

 
};

#endif /* GLOBALWINDOW_H */
