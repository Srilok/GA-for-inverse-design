
def obj_slope(target,predicted,k=5,s=500,dlevel=1):

    import sys
    sys.path.append('/home/share/cnm50256/WTworkflow/scripts/henry')

    import ho,filters
    rec = ho.setup()
#    rec.reportfile = 'report.txt'


    # 'Nsign'     : number of values with wrong signs
    # 'Nwayoff'   : number of values that are wayoff
    #                                                                                                                                                                                                                                   
    # 'D'         : absolute difference
    #               between predicted and target
    # 'MD'        : maximum 'D'
    # 'WMD'       : 'MD' of wayoff values
    # 'AD'        : average relative absolute difference 
    #               between predicted and target
    # 'WAD'       : 'AD' of wayoff values
    # 'SD'        : average relative standard deviation
    #               of predicted values

    #cutoff={'wayoff':0.1}
    #limit = {'Nsign':100,'Nwayoff':100}
    #if filters.slope_slope2(rec,target,predicted,cutoff=cutoff,limit=limit,k=k,s=s,dlevel=dlevel):
    if filters.slope_slope2(rec,target,predicted,k=k,s=s,dlevel=dlevel):
        return rec.objfinal
    return rec.obj
