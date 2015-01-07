ALL_PROGRAMS_SOURCE = \"""Internal;Short name;Full program name;Saturating;Par. 1 name;Par. 1 pop;Par. 2 name;Par. 2 pop;Par. 3 name;Par. 3 pop;Parameter 1;Parameter 2;Parameter 3;COND;Condoms;Condom promotion and distribution;1;Condom usage probability, regular partnerships;All;Condom usage probability, casual partnerships;All;;;M.condom.reg[:];M.condom.cas[:];;SBCC;SBCC;Social and behavior change communication;1;Condom usage probability, regular partnerships;All;Condom usage probability, casual partnerships;All;;;M.condom.reg[:];M.condom.cas[:];;STI;STI;Diagnosis and treatment of sexually transmitted infections;1;Discharging STI prevalence;All;Ulcerative STI prevalence;All;;;M.stiprevdis[:];M.stiprevulc[:];;VMMC;VMMC;Voluntary medical male circumcision;0;Number of circumcisions performed per year;All;;;;;M.numcircum[:];;;FSWP;FSW programs;Programs for female sex workers and clients;1;Condom usage probability, commercial partnerships;FSW;Condom usage probability, commercial partnerships;CSW;HIV testing rates;FSW;M.condom.com[FSW,CSW];M.hivtest[FSW];;MSMP;MSM programs;Programs for men who have sex with men;1;Condom usage probability, regular partnerships;MSM;Condom usage probability, casual partnerships;MSM;;;M.condom.reg[MSM];M.condom.cas[MSM];;PWIDP;PWID programs;Programs for people who inject drugs;1;HIV testing rates;PWID;Condom usage probability, regular partnerships;PWID;Condom usage probability, casual partnerships;PWID;M.hivtest[PWID];M.condom.reg[PWID];M.condom.cas[PWID];OST;OST;Opiate substitution therapy;0;Number of people on OST;All;;;;;M.numost;;;NSP;NSP;Needle-syringe program;1;Needle-syringe sharing rate;All;;;;;M.sharing;;;PREP;PrEP;Pre-exposure prophylaxis;1;PrEP prevalence;All;;;;;M.prep;;;PEP;PEP;Post-exposure prophylaxis;1;PEP prevalence;All;;;;;M.pep;;;CT;Cash transfers;Cash transfers for HIV risk reduction;1;Number of acts per person per year, regular;All;Number of acts per person per year, casual;All;;;M.numacts.reg[:];M.numacts.cas[:];;HTC;HTC;HIV testing and counseling;1;HIV testing rates;All;;;;;M.hivtest[:];;;ART;ART;Antiretroviral therapy;0;Number of people on 1st-line treatment;All;Number of people on 2nd-line treatment;All;;;M.tx1;M.tx2;;PMTCT;PMTCT;Prevention of mother-to-child transmission;0;Number of women on PMTCT;All;;;;;M.numpmtct;;;CARE;Other care;Other care;0;;;;;;;;;;OVC;OVC;Orphans and vulnerable children;0;;;;;;;;;;MGMT;MGMT;Management;0;;;;;;;;;;HR;HR;HR and training;0;;;;;;;;;;ENV;ENV;Enabling environment;0;;;;;;;;;;SP;SP;Social protection;0;;;;;;;;;;MESR;M&E;Monitoring, evaluation, surveillance, and research;0;;;;;;;;;;INFR;INFR;Health infrastructure;0;;;;;;;;;;"""fields = {1:"short_name",2:"name", 3:"saturating"}param_cols = [10,11,12]def programs():    import re    from parameters import maybe_bool    result = []    lines = [l.strip() for l in ALL_PROGRAMS_SOURCE.split('\n')][2:-1]    split_lines = [l.split(';') for l in lines]    for line in split_lines:        entry = dict([(fields[key], maybe_bool(line[key]) ) for key in fields])        params = []        for col in param_cols:            if line[col]:                param, pops = re.match('M\.([^[]+)(?:\[(.+?)\])?',line[col]).groups()                if pops is None: pops=''                params.append({'signature': param.split('.'), 'pops':pops.replace(':','').split(',')})        entry['parameters'] = params        result.append(entry)    return resultprogram_list = programs()def programs_for_input_key(key): #params is the output of parameters.parameters() method    from parameters import input_parameter    param = input_parameter(key)    result = set()    keys = None    if param is not None: keys = param.get('keys')    if keys is not None:        for program in program_list:            for parameter in program['parameters']:                if parameter['signature']==keys:                    result.add(program['short_name'])                    break    return result