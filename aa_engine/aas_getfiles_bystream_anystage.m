function fname = aas_getfiles_bystream_anystage(aap,stage,domain,indices,stream)
[name, ind] = strtok_ptrn(stage,'_0');
index = sscanf(ind,'_%d');
stages = {aap.tasklist.main.module.name};
istage = find(strcmp(stages,name));
aap = aas_setcurrenttask(aap,istage(index));
fname = aas_getfiles_bystream(aap,domain,indices,stream);