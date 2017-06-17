Search.setIndex({docnames:["csg","csg.genetics","csg.genetics.ld","csg.genome","csg.gwas","csg.gwas.Harmonize","csg.intervaltree","csg.notify","csg.pandas","csg.pedigree","csg.pedigree.trios","csg.plotting","csg.statistics","csglib_modules","csglib_tools","index","modules"],envversion:52,filenames:["csg.rst","csg.genetics.rst","csg.genetics.ld.rst","csg.genome.rst","csg.gwas.rst","csg.gwas.Harmonize.rst","csg.intervaltree.rst","csg.notify.rst","csg.pandas.rst","csg.pedigree.rst","csg.pedigree.trios.rst","csg.plotting.rst","csg.statistics.rst","csglib_modules.rst","csglib_tools.rst","index.rst","modules.rst"],objects:{"":{csg:[0,0,0,"-"]},"csg.genetics":{alleles:[1,0,0,"-"],hwe:[1,0,0,"-"],ld:[2,0,0,"-"]},"csg.genetics.alleles":{flip_effects:[1,1,1,""],flip_strand:[1,1,1,""]},"csg.genetics.hwe":{hwe_exact_test:[1,1,1,""]},"csg.genetics.ld":{pyld:[2,0,0,"-"],pyld_example:[2,0,0,"-"]},"csg.genetics.ld.pyld":{Haplotypes:[2,2,1,""],LD:[2,2,1,""]},"csg.genetics.ld.pyld.Haplotypes":{alt:[2,3,1,""],chrom:[2,3,1,""],haplotypes:[2,3,1,""],position:[2,3,1,""],ref:[2,3,1,""],size:[2,3,1,""]},"csg.genetics.ld.pyld.LD":{add_vcf:[2,4,1,""],chunk_M:[2,3,1,""],compute_r:[2,4,1,""],compute_r_matrix:[2,4,1,""],compute_region_freq:[2,4,1,""],compute_variant_freq:[2,4,1,""],get_region_haplotypes:[2,4,1,""],get_variant_haplotypes:[2,4,1,""],get_variant_haplotypes_strict:[2,4,1,""],incr_M:[2,3,1,""],release_vcfs:[2,4,1,""]},"csg.genome":{coordinates:[3,0,0,"-"]},"csg.genome.coordinates":{assign_locus:[3,1,1,""],chrom_int:[3,1,1,""],parse_bp:[3,1,1,""]},"csg.gwas":{Harmonize:[5,0,0,"-"],gwas_catalog:[4,0,0,"-"],independent:[4,0,0,"-"],inflation:[4,0,0,"-"]},"csg.gwas.Harmonize":{Custom:[5,0,0,"-"],EPACTS:[5,0,0,"-"],ImputationQuality:[5,0,0,"-"],Panel:[5,0,0,"-"],QuickTest:[5,0,0,"-"],SNPTEST:[5,0,0,"-"]},"csg.gwas.Harmonize.Custom":{duplicates:[5,1,1,""],harmonize:[5,1,1,""]},"csg.gwas.Harmonize.EPACTS":{duplicates:[5,1,1,""],harmonize_emmax:[5,1,1,""],harmonize_firth:[5,1,1,""]},"csg.gwas.Harmonize.ImputationQuality":{get_imputation_quality:[5,1,1,""],snpstat_entries_it:[5,1,1,""],vcf_entries_it:[5,1,1,""]},"csg.gwas.Harmonize.Panel":{Alleles:[5,2,1,""],Panel:[5,2,1,""]},"csg.gwas.Harmonize.Panel.Alleles":{alt:[5,3,1,""],ref:[5,3,1,""]},"csg.gwas.Harmonize.Panel.Panel":{get_alleles:[5,4,1,""]},"csg.gwas.Harmonize.QuickTest":{duplicates:[5,1,1,""],harmonize_default:[5,1,1,""]},"csg.gwas.Harmonize.SNPTEST":{duplicates:[5,1,1,""],harmonize_frequentist1:[5,1,1,""]},"csg.gwas.gwas_catalog":{PublishedSNVAssociation:[4,2,1,""],get_known_associations_owl:[4,1,1,""],get_known_haplotype_associations_tsv:[4,1,1,""],get_known_snv_associations_tsv:[4,1,1,""],get_study_by_ancestry:[4,1,1,""],map2dbsnp:[4,1,1,""],map2ucsc:[4,1,1,""]},"csg.gwas.gwas_catalog.PublishedSNVAssociation":{chrom:[4,3,1,""],position:[4,3,1,""],pubmed:[4,3,1,""],pvalue:[4,3,1,""],rsid:[4,3,1,""],trait:[4,3,1,""]},"csg.gwas.independent":{Association:[4,2,1,""],get_independent_associations:[4,1,1,""]},"csg.gwas.independent.Association":{alt:[4,3,1,""],chrom:[4,3,1,""],position:[4,3,1,""],pvalue:[4,3,1,""],ref:[4,3,1,""]},"csg.gwas.inflation":{AssociationStats:[4,2,1,""],inflation_from_effect_se:[4,1,1,""],inflation_from_pvalue:[4,1,1,""]},"csg.gwas.inflation.AssociationStats":{effect:[4,3,1,""],pvalue:[4,3,1,""],se:[4,3,1,""]},"csg.intervaltree":{IntervalTree:[6,0,0,"-"]},"csg.intervaltree.IntervalTree":{IntervalTree:[6,2,1,""],IntervalTreeNode:[6,2,1,""]},"csg.intervaltree.IntervalTree.IntervalTree":{add:[6,4,1,""],ascending:[6,4,1,""],complementary:[6,4,1,""],descending:[6,4,1,""],get_height:[6,4,1,""],get_intervals_count:[6,4,1,""],get_values_count:[6,4,1,""],interval_overlap:[6,4,1,""],k_first:[6,4,1,""],k_last:[6,4,1,""],k_nearest_left:[6,4,1,""],k_nearest_right:[6,4,1,""],merge:[6,4,1,""],nearest_left:[6,4,1,""],nearest_right:[6,4,1,""],point_intersect:[6,4,1,""]},"csg.intervaltree.IntervalTree.IntervalTreeNode":{compare:[6,4,1,""],end:[6,3,1,""],get_grandparent:[6,4,1,""],get_height:[6,4,1,""],get_intervals_count:[6,4,1,""],get_sibling:[6,4,1,""],get_uncle:[6,4,1,""],get_values_count:[6,4,1,""],start:[6,3,1,""],values:[6,3,1,""]},"csg.jupyter":{df_to_uri:[0,1,1,""],in_notebook:[0,1,1,""],show_table:[0,1,1,""],text:[0,1,1,""]},"csg.patsy":{get_variables:[0,1,1,""]},"csg.pedigree":{trios:[10,0,0,"-"]},"csg.pedigree.trios":{trios:[10,0,0,"-"]},"csg.pedigree.trios.trios":{KinshipEntry:[10,2,1,""],Relation:[10,2,1,""],SexEntry:[10,2,1,""],Trio:[10,2,1,""],get_trios:[10,1,1,""],king_reader:[10,1,1,""],kinship2degree:[10,1,1,""],pcrelate_reader:[10,1,1,""],sex_reader:[10,1,1,""],write_trios:[10,1,1,""]},"csg.pedigree.trios.trios.KinshipEntry":{id1:[10,3,1,""],id2:[10,3,1,""],kinship:[10,3,1,""],prob_ibd0:[10,3,1,""],prop_ibs0:[10,3,1,""]},"csg.pedigree.trios.trios.Relation":{degree:[10,3,1,""],kinship:[10,3,1,""],prob_ibd0:[10,3,1,""],prop_ibs0:[10,3,1,""]},"csg.pedigree.trios.trios.SexEntry":{id:[10,3,1,""],sex:[10,3,1,""]},"csg.pedigree.trios.trios.Trio":{father:[10,3,1,""],individual:[10,3,1,""],mother:[10,3,1,""],sex:[10,3,1,""]},"csg.plotting":{matplotext:[11,0,0,"-"],qqplot:[11,0,0,"-"],util:[11,0,0,"-"]},"csg.plotting.matplotext":{add_identity:[11,1,1,""]},"csg.plotting.qqplot":{qqplot:[11,1,1,""]},"csg.plotting.util":{thin_dframe:[11,1,1,""]},"csg.statistics":{FDR:[12,0,0,"-"]},"csg.statistics.FDR":{padjust_bh:[12,1,1,""],padjust_qvalue:[12,1,1,""]},"csg.util":{timer:[0,1,1,""]},"csg.yaml":{read_yaml:[0,1,1,""]},csg:{genetics:[1,0,0,"-"],genome:[3,0,0,"-"],gwas:[4,0,0,"-"],intervaltree:[6,0,0,"-"],jupyter:[0,0,0,"-"],notify:[7,0,0,"-"],pandas:[8,0,0,"-"],patsy:[0,0,0,"-"],pedigree:[9,0,0,"-"],plotting:[11,0,0,"-"],statistics:[12,0,0,"-"],util:[0,0,0,"-"],yaml:[0,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","function","Python function"],"2":["py","class","Python class"],"3":["py","attribute","Python attribute"],"4":["py","method","Python method"]},objtypes:{"0":"py:module","1":"py:function","2":"py:class","3":"py:attribute","4":"py:method"},terms:{"1000g":13,"case":5,"class":[2,4,5,6,10],"float":1,"function":[0,1],"int":[3,6,10],"long":6,"new":[0,3,6],"return":[0,1,3,6,10,11],"true":11,Are:11,The:[3,6,11,13,15],There:0,add:6,add_ident:11,add_vcf:2,added:3,addit:13,algorithm:0,alia:[2,4,5,10],all:[0,1,6,13],allel:[0,5,13,14,16],allelea_field:5,alleleb_field:5,alreadi:11,already_transform:11,als:[],also:15,alt:[2,4,5],alt_idx:5,analys:15,ancestri:4,ani:[6,13],api:13,arg:0,arrai:[1,11],ascend:[6,13],assign:3,assign_locu:3,associ:[4,6,15],association_stat:4,associationstat:4,assum:1,attempt:0,axes:11,base:[2,3,4,5,6,10,13],befor:11,between:[3,13],bgzip:13,binari:[6,13],black:[6,13],both:6,branch:6,browser:0,call:[0,3],can:[3,11],cannot:1,catalog_ancestry_fil:4,catalog_associations_fil:4,check:6,cherri:15,chr9:3,chrom:[2,3,4,5],chrom_field:5,chrom_idx:5,chrom_int:3,chrome:0,chromosom:[3,6,13],chromx:3,chrx:3,chunk_m:2,click:0,closest:6,code:[0,15],coded_allele_field:5,coded_allele_freq_field:5,col_chrom:3,col_po:3,color:6,column:[3,11],com:[0,1],comment:5,compar:6,complementari:[6,13],complex:[6,13],compress:13,comput:13,compute_r:2,compute_r_matrix:2,compute_region_freq:2,compute_variant_freq:2,construct:[6,13],contain:3,content:[14,16],control:5,convers:3,convert:[0,3],coordin:[0,14,16],copi:15,cover:6,creat:[0,11],csg:[14,15],csglib_modul:[],current:[0,6],custom:[0,4],cyallel:[0,14,16],data:[0,3,6,11,13],datafram:3,datastuctur:6,dbsnp_file:4,degre:10,deriv:10,descend:[6,13],detail:13,determin:0,deviat:[],df_to_uri:0,dframe:[0,11],dictionari:11,digit:11,discov:0,disequilibrium:13,displai:0,dist:3,distanc:[3,6],document:13,doe:0,done:6,down:11,download:0,duplic:5,each:3,effeci:6,effect:[1,4],effect_field:5,effici:[6,13],either:3,element:0,encod:0,end:[3,6],end_posit:2,entri:5,epact:[0,4],eqtl:0,even:0,everi:6,exact:1,exampl:[0,3],except:10,execut:0,expect:11,extract:[0,13],fals:11,father:10,fdr:[0,14,16],field:[2,4,5,10],file:13,filenam:0,filepath:0,find:[6,11,13],first:13,flip:1,flip_effect:1,flip_strand:1,follow:13,format:[0,13],formula:0,frame:[0,3,11],frequenc:13,from:[3,6,10,13,15],fstr:3,further:13,gap:[6,13],genet:[0,10,14,16],genom:[0,14,15,16],get:13,get_allel:5,get_grandpar:6,get_height:6,get_imputation_qu:5,get_independent_associ:4,get_intervals_count:6,get_known_associations_owl:4,get_known_haplotype_associations_tsv:4,get_known_snv_associations_tsv:4,get_region_haplotyp:2,get_sibl:6,get_study_by_ancestri:4,get_trio:10,get_uncl:6,get_values_count:6,get_vari:0,get_variant_haplotyp:2,get_variant_haplotypes_strict:2,given:[0,1,11,13],grandpar:6,graph:4,greater:10,group:15,gwa:[0,14,15,16],gwas_catalog:[0,14,16],gwascatalog:15,haplotyp:2,haplotypes1:2,haplotypes2:2,hardy_weinberg_equilibrium_exact_test:1,harmon:[0,4,15],harmonize_default:5,harmonize_emmax:5,harmonize_firth:5,harmonize_frequentist1:5,has:6,hdf:[0,14,16],height:6,helper:0,html:0,http:[0,1],hwe:[0,14,16],hwe_exact_test:1,id1:10,id2:10,ident:6,implement:[6,13],imputation_fil:5,imputationqu:[0,4],in_fil:5,in_kinship_fil:10,in_notebook:0,in_sex_fil:10,in_vcf:4,includ:[6,15],incr_m:2,independ:[0,14,16],index:[1,13],individu:10,inflat:[0,14,16],inflation_from_effect_s:4,inflation_from_pvalu:4,inflationfromfil:15,inlin:0,insert:3,integ:3,intersect:6,interv:[0,6,13,14,16],interval_overlap:6,interval_tre:[],intervaltre:[0,14,15,16],intervaltreenod:6,intervaltreenond:[],intev:6,invnorm:[0,14,16],ipython:0,iter:6,iterhelp:[14,16],its:6,jupyt:[14,16],just:[3,13],k_first:6,k_last:6,k_nearest_left:6,k_nearest_right:6,king_read:10,kinship2degre:10,kinship:10,kinship_entri:10,kinshipentri:10,kwarg:0,kwd:0,label:11,larg:0,last:13,left:[6,13],left_rang:5,leftmost:[6,13],length:3,level:6,line_arg:11,line_kwarg:11,link:0,link_text:0,linkag:13,list:[6,11],locu:3,locus_nam:3,log10:11,log:[6,13],longest:6,mai:15,map2dbsnp:4,map2ucsc:4,match:1,math:[14,16],matplotext:[0,14,16],max_end:6,max_locus_bp:4,max_s:5,maxim:6,measur:6,merg:[6,13],min_degre:10,min_info:5,min_mac:5,modifi:1,modul:[14,16],mother:10,na_rep:0,name:[0,3,11],nan:1,ndarrai:1,nearest:13,nearest_left:6,nearest_right:6,need:13,neqtl:0,node:6,non:13,noncoded_allele_field:5,none:[4,6,10],notebook:0,noth:1,notifi:[0,14,16],number:[2,3,4,5,6,10,11],numpi:11,obs_het:1,obs_hom1:1,obs_hom2:1,one:6,ones:0,onli:[0,11],only_col:11,oper:[0,13],option:6,order:[6,13],origin:11,other:[1,15],out:0,out_fil:[5,10],over:6,overlap:[6,13],packag:[14,15,16],padjust_bh:12,padjust_qvalu:12,page:[],pair:[3,13],pairwis:13,panda:[0,11,14,16],panel:[0,4],panel_vcf:5,paramet:[1,3,6,11],pars:3,parse_bp:3,patsi:[14,16],pcrelate_read:10,pedigre:[0,14,15,16],person:15,phase:13,php:1,pick:15,pipelin:15,place:1,pleas:13,plot:[0,14,16],point:[6,11,13],point_intersect:6,pos:3,pos_field:5,posit:[2,3,4,5,6],position_field:5,position_idx:5,possibl:3,print:0,prob_ibd0:10,project:15,prop_ibs0:10,provid:13,publishedsnvassoci:4,pubm:4,pushbullet:[0,14,16],pvalu:[4,12],pvalue_field:5,pvalue_threshold:4,pyld:[0,1],pyld_exampl:[0,1],pyld_example_config:[0,1],pypedia:1,python:15,qqplot:[0,14,16],quality_idx:5,quantil:11,queri:[6,13],question:0,quicktest:[0,4],rais:10,read_yaml:0,red:[6,13],ref:[2,4,5],ref_idx:5,refer:13,referenc:6,region:13,rel:1,relat:10,release_vcf:2,repres:6,requir:3,respect:6,reus:[],right:[6,13],right_rang:5,rightmost:[6,13],round:11,routin:13,row:[3,11],rsid:4,rsq_threshold:4,run_algorithm:0,scale:11,se_field:5,search:[],second:6,see:0,sep:[0,5],seri:[1,11],set:1,sex:10,sex_entri:10,sex_read:10,sexentri:10,should:[3,11],show_tabl:0,shown:0,sibl:6,singl:6,size:[1,2],smaller:6,snap:0,snippet:15,snp:13,snpstat_entries_it:5,snptest:[0,4],special:13,specif:15,specifi:[1,3,6,13],stackoverflow:0,standalon:15,start:6,start_posit:2,statist:[0,14,16],store:6,str:3,string:[0,3],structur:[6,13],studi:15,submodul:[9,14,16],subpackag:[14,16],subset:11,support:13,tab:0,tabix:13,tabl:0,target:1,test:1,text:0,than:10,thi:[0,1,6,13,15],thin:11,thin_dfram:11,those:11,time:[0,6,13],timer:0,too:0,tool:[],total:6,toward:1,trait:4,transform:11,travers:13,tree:[6,13],trio:[0,9],trio_entri:10,triopul:15,tupl:[2,4,5,10],turn:6,two:6,type:[3,6,10],ucsc_snp_fil:4,uncl:6,uniqu:11,unit:3,uri:0,use:15,used:15,user:[0,15],using:[0,11,13],utf:0,util:[14,16],valu:[6,10,11],variabl:0,variant:[0,3],variantid_idx:5,variou:15,vcf:[5,13],vcf_entries_it:5,vcf_path:2,veri:[0,15],well:15,were:0,what:3,when:0,which:11,wide:15,within:[0,6],work:0,write_trio:10,www:1,yaml:[14,16],yield:6},titles:["csg package","csg.genetics package","csg.genetics.ld package","csg.genome package","csg.gwas package","csg.gwas.Harmonize package","csg.intervaltree package","csg.notify package","csg.pandas package","csg.pedigree package","csg.pedigree.trios package","csg.plotting package","csg.statistics package","csglib\u2019s modules","csglib\u2019s tools","Welcome to csglib","csg"],titleterms:{allel:1,content:[0,1,2,3,4,5,6,7,8,9,10,11,12],coordin:3,csg:[0,1,2,3,4,5,6,7,8,9,10,11,12,16],csglib:[13,14,15],csglib_tool:[],custom:5,cyallel:1,data:[],document:[],epact:5,fdr:12,genet:[1,2],genom:3,gwa:[4,5],gwas_catalog:4,harmon:5,hdf:8,hwe:1,imputationqu:5,independ:4,indic:[],inflat:4,interv:3,intervaltre:[6,13],invnorm:12,iterhelp:0,jupyt:0,math:0,matplotext:11,modul:[0,1,2,3,4,5,6,7,8,9,10,11,12,13,15],notifi:7,packag:[0,1,2,3,4,5,6,7,8,9,10,11,12],panda:8,panel:5,patsi:0,pedigre:[9,10],plot:11,pushbullet:7,pyld:2,pyld_exampl:2,pyld_example_config:2,qqplot:11,quicktest:5,snptest:5,statist:12,structur:[],submodul:[0,1,2,3,4,5,6,7,8,10,11,12],subpackag:[0,1,4,9],tabl:[],tool:[14,15],tree:[],trio:10,util:[0,11],welcom:15,yaml:0}})