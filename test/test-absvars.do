// -------------------------------------------------------------------------------------------------
// Test a list of valid absvars
// -------------------------------------------------------------------------------------------------

sysuse auto, clear
set more off
cls


local absvar1		turn
local G_1			1
local ivars1_1		turn
local cvars1_1		
local inter1_1		1
local numsl1_1		0

local absvar2		i.turn
local G_2			1
local ivars1_2		turn
local cvars1_2		
local inter1_2		1
local numsl1_2		0

local absvar3		tur
local G_3			1
local ivars1_3		turn
local cvars1_3		
local inter1_3		1
local numsl1_3		0

local absvar4		i.turn trunk
local G_4			2
local ivars1_4		turn
local cvars1_4		
local inter1_4		1
local numsl1_4		0

local ivars2_4		trunk
local cvars2_4		
local inter2_4		1
local numsl2_4		0

local absvar5		trunk i.turn
local G_5			2
local ivars1_5		trunk
local cvars1_5		
local inter1_5		1
local numsl1_5		0

local ivars2_5		turn
local cvars2_5		
local inter2_5		1
local numsl2_5		0		

local absvar6		turn#trunk
local G_6			1
local ivars1_6		turn trunk
local cvars1_6		
local inter1_6		1
local numsl1_6		0

local absvar7		i.turn#trunk
local G_7			1
local ivars1_7		turn trunk
local cvars1_7		
local inter1_7		1
local numsl1_7		0		

local absvar8		turn#i.trunk
local G_8			1
local ivars1_8		turn trunk
local cvars1_8		
local inter1_8		1
local numsl1_8		0

local absvar9		i.turn#i.trunk
local G_9			1
local ivars1_9		turn trunk
local cvars1_9		
local inter1_9		1
local numsl1_9		0	

local absvar10		i.turn#foreign#trunk
local G_10			1
local ivars1_10		turn foreign trunk
local cvars1_10		
local inter1_10		1
local numsl1_10		0

local absvar11		turn for trunk disp mpg
local G_11			5
local ivars1_11		turn
local cvars1_11		
local inter1_11		1
local numsl1_11		0

local ivars2_11		foreign
local cvars2_11		
local inter2_11		1
local numsl2_11		0

local ivars3_11		trunk
local cvars3_11		
local inter3_11		1
local numsl3_11		0

local ivars4_11		displacement
local cvars4_11		
local inter4_11		1
local numsl4_11		0

local ivars5_11		mpg
local cvars5_11		
local inter5_11		1
local numsl5_11		0

local absvar12		turn trunk#c.gear
local G_12			2
local ivars1_12		turn
local cvars1_12		
local inter1_12		1
local numsl1_12		0

local ivars2_12		trunk
local cvars2_12		gear_ratio
local inter2_12		0
local numsl2_12		1

local absvar13		turn trunk##c.gear
local G_13			2
local ivars1_13		turn
local cvars1_13		
local inter1_13		1
local numsl1_13		0

local ivars2_13		trunk
local cvars2_13		gear_ratio
local inter2_13		1
local numsl2_13		1

local absvar14		turn foreign trunk#c.(gear)
local G_14			3

local ivars1_14		turn
local cvars1_14		
local inter1_14		1
local numsl1_14		0

local ivars2_14		foreign
local cvars2_14		
local inter2_14		1
local numsl2_14		0

local ivars3_14		trunk
local cvars3_14		gear_ratio
local inter3_14		0
local numsl3_14		1	

local absvar15		turn##c.(gear weight)
local G_15			1
local ivars1_15		turn
local cvars1_15		gear_ratio weight
local inter1_15		1
local numsl1_15		2

local absvar16		turn trunk#c.(gear weight length)
local G_16			2

local ivars1_16		turn
local cvars1_16		
local inter1_16		1
local numsl1_16		0

local ivars2_16		trunk
local cvars2_16		gear_ratio weight length
local inter2_16		0
local numsl2_16		3

*local absvar17		turn trunk#c.(gear weight length) foreign
local G_17			1
local ivars1_17		
local cvars1_17		
local inter1_17		
local numsl1_17		

local absvar18		turn c.(gear weight length)#trunk
local G_18			2

local ivars1_18		turn
local cvars1_18		
local inter1_18		1
local numsl1_18		0

local ivars2_18		trunk
local cvars2_18		gear_ratio weight length
local inter2_18		0
local numsl2_18		3

local absvar19		turn c.(gear weight length)#i.trunk
local G_19			2

local ivars1_19		turn
local cvars1_19		
local inter1_19		1
local numsl1_19		0

local ivars2_19		trunk
local cvars2_19		gear_ratio weight length
local inter2_19		0
local numsl2_19		3

local absvar20		turn (c.gear c.weight c.length)#trunk
local G_20			2

local ivars1_20		turn
local cvars1_20		
local inter1_20		1
local numsl1_20		0

local ivars2_20		trunk
local cvars2_20		gear_ratio weight length
local inter2_20		0
local numsl2_20		3

local absvar21		turn trunk#c.gear foreign##c.(weight length)
local G_21			3
local ivars1_21		turn
local cvars1_21		
local inter1_21		1
local numsl1_21		0

local ivars2_21		trunk
local cvars2_21		gear_ratio
local inter2_21		0
local numsl2_21		1

local ivars3_21		foreign
local cvars3_21		weight length
local inter3_21		1
local numsl3_21		2

local absvar22		FE1=turn foreign FE3=trunk
local G_22			3
local ivars1_22		turn
local cvars1_22		
local inter1_22		1
local numsl1_22		0
local target1_22 	FE1

local ivars2_22		foreign
local cvars2_22		
local inter2_22		1
local numsl2_22		0

local ivars3_22		trunk
local cvars3_22		
local inter3_22		1
local numsl3_22		0
local target3_22 	FE3

local absvar23		turn FE=foreign#c.gear
local G_23			2
local ivars1_23		turn
local cvars1_23		
local inter1_23		1
local numsl1_23		0

local ivars2_23		foreign
local cvars2_23		gear_ratio
local inter2_23		0
local numsl2_23		1
local target2_23	FE

local absvar24		turn FE=foreign#c.(gear length)
local G_24			2
local ivars1_24		turn
local cvars1_24		
local inter1_24		1
local numsl1_24		0

local ivars2_24		foreign
local cvars2_24		gear_ratio length
local inter2_24		0
local numsl2_24		2
local target2_24	FE

local absvar25		turn foreign, savefe
local G_25			2
local savefe_25 	1

local ivars1_25		turn
local cvars1_25		
local inter1_25		1
local numsl1_25		0

local ivars2_25		foreign
local cvars2_25		
local inter2_25		1
local numsl2_25		0


local absvar26		FE1 = turn foreign FE3 = trunk
local G_26			3
local ivars1_26		turn
local cvars1_26		
local inter1_26		1
local numsl1_26		0
local target1_26 	FE1

local ivars2_26		foreign
local cvars2_26		
local inter2_26		1
local numsl2_26		0

local ivars3_26		trunk
local cvars3_26		
local inter3_26		1
local numsl3_26		0
local target3_26 	FE3

local absvar27		"FE1 	=	 turn foreign FE3	=trunk"
local G_27			3
local ivars1_27		turn
local cvars1_27		
local inter1_27		1
local numsl1_27		0
local target1_27 	FE1

local ivars2_27		foreign
local cvars2_27		
local inter2_27		1
local numsl2_27		0

local ivars3_27		trunk
local cvars3_27		
local inter3_27		1
local numsl3_27		0
local target3_27 	FE3

// -------------------------------------------------------------------------------------------------

set trace off
forval i = 1/30 {
	local absvar `absvar`i''
	local G_1			1
	if ("`absvar'"=="") continue
	di as input "[`i'] `absvar'"
	ParseAbsvars `absvar'
	return list

	local G = `G_`i''
	assert_msg r(G)==`G', msg("r(G)=`r(G)', expected `G'")
	assert_msg "`r(savefe)'"=="`savefe_`i''" | ("`r(savefe)'"=="0" & "`savefe_`i''"==""), ///
		msg("r(savefe)=`r(savefe)', expected `savefe_`i''")
	forval g=1/`G' {
		assert_msg "`r(target`g')'"=="`target`g'_`i''" , msg(`""`r(target`g')'"=="`target`g'_`i''""')
		assert_msg "`r(ivars`g')'"=="`ivars`g'_`i''" , msg(`""`r(ivars`g')'"=="`ivars`g'_`i''""')
		assert_msg "`r(cvars`g')'"=="`cvars`g'_`i''" , msg(`""`r(cvars`g')'"=="`cvars`g'_`i''""')
		assert_msg `r(has_intercept`g')'==`inter`g'_`i'' , msg(`"`r(has_intercept`g')'==`inter`g'_`i''"')
		assert_msg `r(num_slopes`g')'==`numsl`g'_`i'' , msg(`"`r(num_slopes`g')'==`numsl`g'_`i''"')
	}
}

* TODO: Check that we flag the most common errors

exit
