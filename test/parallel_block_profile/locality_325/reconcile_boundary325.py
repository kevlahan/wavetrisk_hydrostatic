"""Physical closure of the live fixed-Domain to block-boundary import route.

The bulk writeback manifest runs block->Domain; boundary numeric traffic runs
in reverse. Rows are matched to actual installed boundary storage, including
retained local blocks. No rank-label equality is assumed for Domain providers.
"""
def require(ok,code):
    if not ok:raise ValueError('boundary325.'+code)

def prefix(xs):
    out=[];total=0
    for x in xs:out.append(total);total+=x
    return out

def reconcile_boundary325(plans,contexts,owner_of):
    """owner_of receives a one-based physical catalogue ID. Returns invariant rows."""
    n=len(plans);require(n==len(contexts),'rank_extent');pre='block_writeback_plan%'
    blocks={};rows=[]
    for rank,ctx in enumerate(contexts):
        for block in ctx['blocks']:
            cat=block['catalog'];f=block['fields'];ident=f['identity'];raw=f['bdry_storage']
            require(1<=cat<=640 and cat not in blocks and owner_of(cat)==rank,'block_owner_census')
            require(len(ident)==13 and ident[0]==cat-1 and ident[1]==(cat-1)//4,'block_identity')
            require(len(raw)%7==0,'storage_extent')
            boundaries=[raw[i:i+7] for i in range(0,len(raw),7)]
            ns=ident[5]*ident[8]*ident[9];nv=ident[8]*ident[10]
            require(ns>0 and nv>0,'field_width')
            for i,b in enumerate(boundaries,1):
                require(b[0]>=0 and b[1]>=0 and b[2]>=0 and b[5]>0 and b[6]>=0,'storage_bounds')
                rows.append((cat,i,(cat-1)//4,b[0],b[1],b[2],b[5],b[6],ns*b[5],nv*b[5]))
            blocks[cat]=(rank,boundaries,ns,nv)
    require(set(blocks)==set(range(1,641)),'complete_block_census')
    for rank,p in enumerate(plans):
        manifests={}
        for direction in ('send','recv'):
            num=p['n_'+direction][0]
            cats=p[pre+direction+'_block'][:num]
            expected=sorted([cat for cat in range(1,641)
                      if (owner_of(cat)==rank and (cat-1)//4!=rank if direction=='send'
                          else (cat-1)//4==rank and owner_of(cat)!=rank)],
                      key=lambda cat:((cat-1)//4 if direction=='send' else owner_of(cat),cat))
            require(cats==expected,'bulk_'+direction+'_catalogue')
            counts=[];scalar=[];vector=[];sources=[];levels=[];nodes=[]
            for cat in cats:
                _,bb,ns,nv=blocks[cat]
                counts.append(len(bb));scalar.append(sum(b[5]*ns for b in bb));vector.append(sum(b[5]*nv for b in bb))
                sources.extend(b[0] for b in bb);levels.extend(b[1] for b in bb);nodes.extend(b[5] for b in bb)
            def equal(field,want,empty_fill=0):
                # Source allocates max(1,logical_size); unused singleton is initialized.
                actual=p[pre+direction+'_boundary_'+field]
                require(actual==(want if want else [empty_fill]),direction+'_'+field)
            equal('count',counts);equal('displ',prefix(counts))
            equal('scalar_nvalue',scalar);equal('vector_nvalue',vector)
            equal('source',sources);equal('level',levels,-1);equal('nnode',nodes)
            require(p['n_'+direction+'_boundary']==[len(sources)],direction+'_total')
            manifests[direction]=(cats,scalar,vector)
        for component,column in (('scalar',1),('vector',2)):
            for flow,direction in (('domain_send','recv'),('block_recv','send')):
                cats=manifests[direction][0];values=manifests[direction][column]
                count=[0]*n
                for cat,value in zip(cats,values):
                    peer=owner_of(cat) if flow=='domain_send' else (cat-1)//4
                    count[peer]+=value
                field=pre+'boundary_'+component+'_'+flow
                require(p[field+'_count']==count,component+'_'+flow+'_count')
                require(p[field+'_displ']==prefix(count),component+'_'+flow+'_displ')
                require(p[field+'_buffer.shape']==[max(1,sum(count))],component+'_'+flow+'_capacity')
        retained=sum(owner_of(cat)==rank==(cat-1)//4 for cat in range(1,641))
        require(p['n_retained']==[retained],'retained_coverage')
    return sorted(rows)
