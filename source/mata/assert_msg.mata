mata:
mata set matastrict on
	void assert_msg(real scalar t, | string scalar msg)
	{
		if (args()<2 | msg=="") msg = "assertion is false"
	        if (t==0) _error(msg)
	}
end
