// Used inside the connected presentation session. Evidence comes from the saved R preparation script.
export function buildSummarySlide(presentation, e) {
  if(e.policy!=='observed_only' || e.volcanoes.length!==3) throw new Error('Expected three observed-only interactions.');
  for(const v of e.volcanoes) if(v.series.reduce((n,s)=>n+s.xvalues.length,0)!==1858) throw new Error('Wrong volcano universe.');
  const ids=presentation.slides.items.map(s=>s.id);
  if(ids.length) presentation.slides.delete(ids);
  const slide=presentation.slides.add({width:1600,height:900});
  slide.background.fill='#FFFFFF';
  const ink='#102C45', muted='#526577', orange='#D55E00', blue='#0072B2', light='#C7CFD6';
  const excelNumber=n=>Number(n.toPrecision(15)); // Excel-compatible chart display precision; source statistics stay unchanged.
  function text(name,value,x,y,w,h,size=24,bold=false,color=ink) {
    const shape=slide.shapes.add({name,geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
    shape.text=value;
    shape.text.style={fontFamily:'Arial',fontSize:size,bold,color,verticalAlignment:'middle'};
    return shape;
  }
  text('title','JE2 PASN responses: progress through Step 5',40,24,1520,62,46,true);
  text('workflow','October total proteome: 27 samples across three genotypes and three media.\nUniProt mapping retained 1,858 source-observed proteins for interaction testing and GO/KEGG enrichment.',40,91,1520,68,23,false,muted);
  text('volcano-heading','Protein abundance interactions',40,166,1520,34,27,true);
  const fills={other:light,negative:blue,positive:orange};
  const axisText={typeface:'Arial',fontSize:17,fill:muted};
  e.volcanoes.forEach((v,i)=>{
    const series=v.series.map(s=>({name:s.name,xValues:s.xvalues.map(excelNumber),values:s.yvalues.map(excelNumber),
      fill:fills[s.name],line:{fill:'none',width:0},marker:{symbol:'circle',size:s.name==='other'?2:4,fill:fills[s.name]}}));
    for(const x of [-1,1]) series.push({name:'Effect cutoff '+x,xValues:[x,x],values:[0,e.y_limit],
      line:{fill:'#AEB9C4',width:.7,style:'dashed'},marker:{symbol:'none'}});
    series.push({name:'FDR cutoff',xValues:[-e.x_limit,e.x_limit],values:[1,1],
      line:{fill:'#AEB9C4',width:.7,style:'dashed'},marker:{symbol:'none'}});
    const c=slide.charts.add('scatter',{
      position:{left:35+i*520,top:205,width:505,height:278},title:v.title,
      titleTextStyle:{typeface:'Arial',fontSize:25,bold:true,fill:ink},titlePlacement:'aboveChart',
      scatterOptions:{style:'lineWithMarkers'},series,hasLegend:false,
      chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',chartLine:{fill:'none',width:0},plotAreaLine:{fill:'none',width:0},
      xAxis:{title:{text:'Interaction effect (log2)',textStyle:axisText},min:-e.x_limit,max:e.x_limit,majorUnit:4,
        textStyle:axisText,position:'bottom',tickLabelPosition:'low',line:{fill:'#AEB9C4',width:1},majorGridlines:null},
      yAxis:{title:{text:'-log10(protein FDR)',textStyle:axisText},min:0,max:e.y_limit,majorUnit:3,
        textStyle:axisText,position:'left',tickLabelPosition:'low',line:{fill:'#AEB9C4',width:1},majorGridlines:null}
    });
    text('hit-counts-'+i,'Selected hits: +'+v.counts.positive+' / -'+v.counts.negative,105+i*520,480,410,25,19,false,muted);
  });
  text('colour-key','Orange: positive interaction. Blue: negative. Hits require protein FDR <0.1 and |effect| >1.',40,510,1520,27,19,false,muted);
  text('gsea-heading','Selected KEGG GSEA modules',40,553,760,37,27,true);
  const labels=['Purine metabolism','Oxidative phosphorylation','Nitrogen metabolism','Aminoacyl-tRNA biosynthesis'];
  const values=[['Signed NES','stk / anc','stp / anc','stp / stk'],...e.gsea.rows.map((r,i)=>[labels[i],...r.cells.map(c=>c.significant?(c.NES>0?'+':'')+c.NES.toFixed(2):'ns')])];
  const table=slide.tables.add({rows:5,columns:4,left:40,top:595,width:760,height:205,
    columnWidths:[364,132,132,132],values});
  table.cells.block({row:0,column:0,rowCount:5,columnCount:4}).assign({margins:{top:3,bottom:3,left:8,right:8},anchor:'center'});
  table.borders.assign({fill:'#FFFFFF',width:2,style:'solid'});
  for(let r=0;r<5;r++) {
    table.rows[r].height=41;
    for(let col=0;col<4;col++) {
      const cell=table.getCell(r,col);
      cell.text.style={fontFamily:'Arial',fontSize:r===0?20:21,bold:r===0||col===0,color:r===0?'#FFFFFF':ink};
      cell.fill=r===0?ink:(col===0?'#F0F4F7':'#F8FAFC');
      if(r>0&&col>0) {
        const val=e.gsea.rows[r-1].cells[col-1];
        cell.fill=val.significant?(val.NES>0?'#F7D8C3':'#CCE2F1'):'#F8FAFC';
        cell.text.style={fontFamily:'Arial',fontSize:23,bold:val.significant,color:val.significant?ink:'#8B99A5'};
      }
    }
  }
  text('gsea-key','Positive: larger SN-TSB response in second genotype. ns: study FDR >=0.1.',40,805,760,30,17,false,muted);
  text('ora-heading','Selected GO BP ORA terms, stk vs ancestor',850,553,710,37,27,true);
  const ora=e.ora.rows;
  const displayOra=[...ora].reverse(); // Native horizontal bars place the first category at the bottom.
  const shortNames={'GO:0019835':'Cytolysis (stk)','GO:0009082':'BCAA biosynthesis (anc)','GO:0006189':'IMP biosynthesis (anc)'};
  slide.charts.add('bar',{
    position:{left:840,top:590,width:720,height:210},
    categories:displayOra.map(r=>shortNames[r.term_id]),
    series:[{name:'Fold enrichment',values:displayOra.map(r=>excelNumber(r.fold_enrichment)),
      fill:orange,line:{fill:'none',width:0},points:displayOra.map((r,idx)=>({idx,fill:r.effect_sign==='positive'?orange:blue}))}],
    barOptions:{direction:'bar',grouping:'clustered',gapWidth:85},hasLegend:false,
    chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',chartLine:{fill:'none',width:0},
    xAxis:{textStyle:{...axisText,fontSize:20},line:{fill:'none',width:0},majorGridlines:null},
    yAxis:{title:{text:'Fold enrichment',textStyle:{...axisText,fontSize:18}},min:0,max:22,majorUnit:10,
      textStyle:axisText,line:{fill:'#AEB9C4',width:1},majorGridlines:null},
    dataLabels:{showValue:true,position:'outEnd',numberFormatCode:'0.0',textStyle:{typeface:'Arial',fontSize:21,bold:true,fill:ink}}
  });
  text('ora-key','Labels show the more-positive SN-TSB response. All study FDR <0.1.',850,805,710,30,17,false,muted);
  text('conclusion','Purine differences are clearest in stk. Nitrogen-related differences appear in stp.',40,835,1520,34,26,true);
  text('caveat','Protein abundance does not establish pathway activity. The stp strain also shows reduced Stk1.',40,871,1520,23,17,false,muted);
  // Metadata records evidence and selection, not authoring/check commentary.
  slide.speakerNotes.text=[
    'Scope: protein-abundance work through Step 5, October export. 27 experimental samples comprise 3 genotypes x3 media x3 replicates. SNpulse belongs to the design but is outside these three focal interactions.',
    'Source-observed-only universe: 1858 unique quality-mapped UniProt proteins per interaction. Some observed proteins have only1/3 measured values. Source protein effects, source FDRs and source statistic are unchanged. No source differential model was refitted.',
    'Volcanoes: effect=(A_TSB-B_TSB)-(A_SN-B_SN)=R_B-R_A, R=SN-TSB. Positive is more positive response in the second genotype. All1858 proteins appear in each chart. Selection: source proteinFDR<0.1 andabs(effect)>1. Threshold lines at +/-1 and-log10(.1)=1. Source FDRzeros='+e.zero_FDR_count+'. Chart coordinates use15 significant digits for Excel-compatible data editing; full precision source statistics remain in the saved results.',
    'KEGG heatmap: four selected modules representing nucleotide metabolism, respiration, nitrogen response and translation. Values are signedNES; ns means not passingstudyFDR<0.1. The completeuncappedplots/tables remainavailable. Selection is for thisone-slide overview, not a new testing family.',
    'GO_BP ORA: three nonredundant examples among eight study-significant stk/ancestor terms: '+ora.map(r=>r.term+';fold='+r.fold_enrichment+';overlap='+r.overlap+';studyFDR='+r.FDR+';'+r.direction).join('\n'),
    'Study BH families: methods separately across ten comparisons and five original collections. Effective set sizes10-500. The fixed courseSTRING-v12 snapshot supplies memberships transferred through singletonUniProt mappings.',
    'Caveats: morepositive can mean a largerincrease or smallerdecrease. Enrichment is not a directmeasure of enzymeactivity, secretion ormetabolicflux. The stp strain alsohasreducedStk1, limiting anisolated phosphatase mechanism. BCAA ORA doesnot imply whole-pathway GSEA agreement.',
    'Evidence: results/step04/total/annotated_all.tsv.gz; results/step05/observed_only/{ora_all,gsea_all}.tsv.gz; reports/step05/observed_only/validation.txt. Reproduce evidence with scripts/05_figures/05_prepare_slide_summary.R; authoring with scripts/06_presentation/05_build_summary_slide.mjs. Full report: reports/step05/observed_only/README.md.'
  ].join('\n\n');
  return slide;
}
