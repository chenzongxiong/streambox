/* basically same thing. need a different class for dispatching eval  */
template <typename InputT>
class WindowsHolder : public WindowedSum<InputT> {
	WindowsHolder(string name) : WindowedSum<InputT>(name) { }
};
