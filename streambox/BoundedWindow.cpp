#include "BoundedWindow.h"
// xzl: they're ptime. must be initialized out of line. and cannot be in
// the headers
const ptime BoundedWindow::TIMESTAMP_MIN_VALUE(min_date_time);
const ptime BoundedWindow::TIMESTAMP_MAX_VALUE(max_date_time);
ptime GlobalWindow::END_OF_GLOBAL_WINDOW = TIMESTAMP_MAX_VALUE - days(1);
//GlobalWindow INSTANCE;
