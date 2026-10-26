"""Check accepted data, receipts and numerical gates with numerical entry points blocked."""
import json
import adapter, run_pilot
from adapter import ROOT,PUBLIC,mp,numeric,read,require,digest
from source import xlsx_rows

def main():
    calls=[]
    def forbidden(*args,**kwargs):
        calls.append(True)
        raise RuntimeError('Numerical work is forbidden in saved-record verification')
    for module in [adapter,run_pilot]:
        module.trajectory=forbidden;module.preparation=forbidden
    for name in ['trajectory','evaluate','scale_values','matrix_state']:
        setattr(numeric,name,forbidden)
    for name in ['qr_solve','lu_solve','eig','eigsy','svd']:setattr(mp,name,forbidden)
    cfg=read(ROOT/'config.json');checked=run_pilot.verify(cfg)
    path=PUBLIC/'data/raw/data_akbilgic.xlsx'
    require(digest(path)==cfg['source_sha256'],'Raw source changed')
    rows=xlsx_rows(path,'orjinal data');header=dict(rows)
    require(header[1][2][0]=='USD BASED','Response currency label')
    require([header[2][k][0] for k in range(3,10)]==['SP','DAX','FTSE','NIKKEI','BOVESPA','EU','EM'],'Feature header')
    observed=[(n,c) for n,c in rows if 3<=n<=538]
    require([n for n,_ in observed]==list(range(3,539)),'Source rows')
    tokens=read(ROOT/'input_tokens.json')
    require(tokens['X']==[[c[k][0] for k in range(3,10)] for _,c in observed],'Feature tokens')
    require(tokens['y']==[c[2][0] for _,c in observed],'Response tokens')
    require(all(c[k][1]=='n' and not c[k][2] for _,c in observed for k in range(2,10)),'Numeric source cells')
    require(not calls,'Unexpected numerical work')
    checked.update(raw_feature_cells_verified=3752,raw_response_cells_verified=536,
        new_optimizer_updates=0,new_QR_solves=0,new_eigensolves=0)
    print(json.dumps(checked,indent=2))
    return checked

if __name__=='__main__':main()
