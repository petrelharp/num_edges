import tskit
import pandas as pd
from remove_isolated_unary import remove_isolated_unary
import tscompare
import tsdate
import msprime
import tsinfer

# Data with varying sample
samplelist = [10, 50, 100, 500, 1000]
sample_arf = pd.DataFrame(columns = samplelist, index=['S', 'SE', 'I', 'IE', 'IS', 'ISE', 'ID', 'IDE', 'IDS', 'IDSE'])
sample_tpr = pd.DataFrame(columns = samplelist, index=['S', 'SE', 'I', 'IE', 'IS', 'ISE', 'ID', 'IDE', 'IDS', 'IDSE'])

for sample in samplelist:
    # generate true tree seqeunce ts and inferred tree sequence infer_ts
    ts = msprime.sim_ancestry(sample, population_size=10000, sequence_length=5e7,
                            recombination_rate=1e-8, coalescing_segments_only=False)
    mu = 1.29e-8  # mutation rate for the inferred tree sequences
    ts = msprime.sim_mutations(ts, mu)
    infer_ts = tsinfer.infer(tsinfer.SampleData.from_tree_sequence(ts))
    t = infer_ts.tables
    t.compute_mutation_times()
    infer_ts = t.tree_sequence()
    infer_dated = tsdate.date(tsdate.util.split_disjoint_nodes(infer_ts), 
                      mutation_rate=mu,
                      allow_unary=True,
                      rescaling_intervals=100
                     )
    iid = remove_isolated_unary(infer_dated)
    ts = remove_isolated_unary(ts)
    ts.dump(f'trees2/ts_{sample}s_5e7')
    iid.dump(f'trees2/infer_dated_ts_{sample}s_5e7')
    
    # ts = tskit.load(f'trees/ts2_{sample}s_5e7_nounary')
    # infer_ts = tskit.load(f'trees/infer_ts2_{sample}s_5e7_nounary')

    # compare
    ID = tscompare.compare(iid, ts)
    print(sample, 'ID', ID)
    sample_arf.loc['ID', sample], sample_tpr.loc['ID', sample] = ID.arf, ID.tpr

    iide = iid.extend_haplotypes()
    IDE = tscompare.compare(iide, ts)
    print(sample, 'IDE', IDE)
    sample_arf.loc['IDE', sample], sample_tpr.loc['IDE', sample] = IDE.arf, IDE.tpr
    
    iids = iid.simplify()
    IDS = tscompare.compare(iids, ts)
    print(sample, 'IDS', IDS)
    sample_arf.loc['IDS', sample], sample_tpr.loc['IDS', sample] = IDS.arf, IDS.tpr

    iidse = iids.extend_haplotypes()
    IDSE = tscompare.compare(iidse, ts)
    print(sample, 'IDSE', IDSE)
    sample_arf.loc['IDSE', sample], sample_tpr.loc['IDSE', sample] = IDSE.arf, IDSE.tpr
    
    I = tscompare.compare(infer_ts,ts)
    print(sample, 'I', I)
    sample_arf.loc['I', sample], sample_tpr.loc['I', sample] = I.arf, I.tpr

    ie = infer_ts.extend_haplotypes()
    IE = tscompare.compare(ie,ts)
    print(sample, 'IE', IE)
    sample_arf.loc['IE', sample], sample_tpr.loc['IE', sample] = IE.arf, IE.tpr
    
    ss = ts.simplify()
    S = tscompare.compare(ss,ts)
    print(sample, 'S', S)
    sample_arf.loc['S', sample], sample_tpr.loc['S', sample] = S.arf, S.tpr
    
    sse = ss.extend_haplotypes()
    SE = tscompare.compare(sse,ts)
    print(sample, 'SE', SE)
    sample_arf.loc['SE', sample],  sample_tpr.loc['SE', sample] = SE.arf, SE.tpr

    iis = infer_ts.simplify()
    IS = tscompare.compare(iis,ts)
    print(sample,'IS', IS)
    sample_arf.loc['IS', sample], sample_tpr.loc['IS', sample] = IS.arf, IS.tpr

    ise = iis.extend_haplotypes()
    ISE = tscompare.compare(ise,ts)
    print(sample,'ISE', ISE)
    sample_arf.loc['ISE', sample], sample_tpr.loc['ISE', sample] = ISE.arf, ISE.tpr
    
sample_arf.to_csv('figure6-arf-over-sample.csv')
sample_tpr.to_csv('figure6-tpr-over-sample.csv')

# Data with varying sequence length
lengthlist = [1e6, 5e6, 1e7, 3e7, 5e7]
names = ['1e6', '5e6', '1e7', '3e7', '5e7']
length_arf = pd.DataFrame(columns = names, index=['S', 'SE', 'I', 'IE', 'IS', 'ISE', 'ID', 'IDE', 'IDS', 'IDSE'])
length_tpr = pd.DataFrame(columns = names, index=['S', 'SE', 'I', 'IE', 'IS', 'ISE', 'ID', 'IDE', 'IDS', 'IDSE'])

for (length, name) in zip(lengthlist, names):
    # generate true tree seqeunce ts and inferred tree sequence infer_ts
    ts = msprime.sim_ancestry(1000, population_size=10000,
                              sequence_length=length,
                              recombination_rate=1e-8,
                              coalescing_segments_only=False
                             )
    mu = 1.29e-8  # mutation rate for the inferred tree sequences
    ts = msprime.sim_mutations(ts, mu)
    infer_ts = tsinfer.infer(tsinfer.SampleData.from_tree_sequence(ts))
    t = infer_ts.tables
    t.compute_mutation_times()
    infer_ts = t.tree_sequence()
    infer_dated = tsdate.date(tsdate.util.split_disjoint_nodes(infer_ts), 
                      mutation_rate=mu,
                      allow_unary=True,
                      rescaling_intervals=100
                     )
    iid = remove_isolated_unary(infer_dated)
    ts = remove_isolated_unary(ts)
    ts.dump(f'trees2/ts_1000s_{name}')
    iid.dump(f'trees2/infer_dated_ts_1000s_{name}')
    # ts = tskit.load(f'trees/ts_1000s_{name}_nounary')
    # infer_ts = tskit.load(f'trees/infer_ts_1000s_{name}_nounary')

    # compare
    I = tscompare.compare(infer_ts,ts)
    print(name, 'I', I)
    length_arf.loc['I', name], length_tpr.loc['I', name] = I.arf, I.tpr
    ie = infer_ts.extend_haplotypes()
    IE = tscompare.compare(ie,ts)
    print(name,'IE', IE)
    length_arf.loc['IE', name], length_tpr.loc['IE', name] = IE.arf, IE.tpr

    ID = tscompare.compare(iid, ts)
    print(name, 'ID', ID)
    length_arf.loc['ID', name], length_tpr.loc['ID', name] = ID.arf, ID.tpr

    iide = iid.extend_haplotypes()
    IDE = tscompare.compare(iide, ts)
    print(name, 'IDE', IDE)
    length_arf.loc['IDE', name], length_tpr.loc['IDE', name] = IDE.arf, IDE.tpr

    iids = iid.simplify()
    IDS = tscompare.compare(iids, ts)
    print(name, 'IDS', IDS)
    length_arf.loc['IDS', name], length_tpr.loc['IDS', name] = IDS.arf, IDS.tpr

    iidse = iids.extend_haplotypes()
    IDSE = tscompare.compare(iidse, ts)
    print(name, 'IDSE', ID)
    length_arf.loc['IDSE', name], length_tpr.loc['IDSE', name] = IDSE.arf, IDSE.tpr
    
    ss = ts.simplify()
    S = tscompare.compare(ss,ts)
    print(name, 'S', S)
    length_arf.loc['S', name], length_tpr.loc['S', name] = S.arf, S.tpr
    
    sse = ss.extend_haplotypes()
    SE = tscompare.compare(sse, ts)
    print(name, 'SE', SE)
    length_arf.loc['SE',name],  length_tpr.loc['SE', name] = SE.arf, SE.tpr

    iis = infer_ts.simplify()
    IS = tscompare.compare(iis,ts)
    print(name, 'IS', IS)
    length_arf.loc['IS', name], length_tpr.loc['IS', name] = IS.arf, IS.tpr
    
    ise = iis.extend_haplotypes()
    ISE = tscompare.compare(ise,ts)
    print(name, 'ISE', ISE)
    length_arf.loc['ISE', name], length_tpr.loc['ISE', name] = ISE.arf, ISE.tpr
    
length_arf.to_csv('figure6-arf-over-length.csv')
length_tpr.to_csv('figure6-tpr-over-length.csv')