import Oviz from "crux";
import { editorConfig, editorRef } from "./editor";
import template from "./template.bvt";

import { ComplexBoxplot, processBoxData } from "oviz-components/complex-boxplot";
import { GridPlot } from "oviz-components/grid-plot";
import { EditText } from "oviz-components/edit-text";
import {register} from "page/visualizers";
import { rankDict, sortByRankKey } from "utils/bio-info";
import { registerEditorConfig } from "utils/editor";
import { generatePath, plotData } from "./data";
import * as d3 from "d3";
import * as ds from "d3-sankey";
import * as TextSize from "crux/dist/utils/text-size";
import { TableSimplePlugin } from "bootstrap-vue";
import { isFunctionOrConstructorTypeNode, isTemplateExpression, nodeModuleNameResolver, textChangeRangeIsUnchanged } from "typescript";
import { tooltip } from "crux/dist/utils";
import { co2servermap,imageMap } from "./commonInfo"
import { Text } from "crux/src/element";
import { max } from "crux/dist/utils/math";
import { replace } from "utils/replace";
import { GeneGome } from "./genetree";
import { CircosContentLine } from "crux/src/element/circos";


//const MODULE_NAME = "analysis_gene_tree";
const MODULE_NAME = "analysis-gene-tree";

function init() {
    if (!window.gon || window.gon.module_name !== MODULE_NAME) return;

    const {visualizer} = Oviz.visualize({
        el: "#canvas",
        root:new GeneGome(),
        renderer:"svg",
        width: 2000,
        height: 2400,
        // template,
        components: { GridPlot, ComplexBoxplot, EditText},
        data: {
            // treeColormap:{},
            showcommonname:false,
            showgenesymbol:false,
            showprotein:false,
            geneBlockcolorfull:false,
            genecolorpattern:true,
            showimage:true,
            fillOpacitynum:null,
            specisColor:false,
            treewidth :500,
            leafSeparator:" ",
            radicalPart:false,
            horizonalPart:true,
            getlen(text,sig,size){
                return TextSize.measuredTextSize(text,size)[sig]
            },
            getSum(arr,i,sig){
                let sum = 0
                arr.forEach((item,index)=> {
                    index < i? sum  = sum + item.length: null
                    sig == "text"? (index == i? sum = sum + item.length/2-4:null):null
                });
                return sum
            },
            getgeneColor(ele){
                let index = this.geneSymbol_list.indexOf(ele)
                const color = Oviz.color.Color.hsl((index%6)*60, 60+Math.floor((index/6))*10, 60+Math.floor((index/6))*10)
                return color.string
            },
        },
        loadData: {
            Gene_tree: {
            // tree:{
                fileKey: "demoTree", 
                type: "newick",  
                multiple: true,
                loaded(data){
                    // this.data.genetree = data[0]
                    // console.log("tree data:",data)
                    this.data.allDatas = data
                    this.data.genecomeTree = ""
                    //这里需要辨别treeall和filltertree
                    let treefilslist = []
                    let oritreefilslist = []
                    let mainTreename = ""


                    this.selectedFiles.demoTree.forEach((item,index) => {
                        let temparr = item.url.split("/")
                        temparr[temparr.length-1].indexOf("tree.newick") != -1? mainTreename = temparr[temparr.length-1]:null
                        treefilslist.push(temparr[temparr.length-1])
                        oritreefilslist.push(temparr[temparr.length-1])
                    });
                    this.data.mainTreename = mainTreename

                    this.data.genetree = data[oritreefilslist.indexOf(this.data.mainTreename)] // 主


                    //构建供选择的列表 替代用户直接选genesymbol
                    treefilslist = treefilslist.filter(d => d.indexOf("filter")!=-1)
                    // console.log("treefilslist:",treefilslist) //filter子树
                    this.data.treefilslist = treefilslist
                    this.data.oritreefilslist = oritreefilslist
                    

                }
            },
            Species_family:{
            // protein2geneall:{
                fileKey:"protein2geneall",
                type:"tsv",
                loaded(data){
                    // console.log("protein2geneall:",data)

                    this.data.orifamilyinfo = data
                    // let geneSymbol_list = Array.from(new Set(data.map(x=>x["Gene_symbol"]))).filter(x=>x!= "").filter(x=>x!=undefined)
                    let geneSymbol_list = Array.from(new Set(data.map(x=>x["Gene"]))).filter(x=>x!= "").filter(x=>x!=undefined)
                    // let proteinid_list = Array.from(new Set(data.map(x=>x["Protein_ID"]))).filter(x=>x!= "").filter(x=>x!=undefined)
                    let proteinid_list = Array.from(new Set(data.map(x=>x["Best_Peptide"]))).filter(x=>x!= "").filter(x=>x!=undefined)
                    // let species_list = Array.from(new Set(data.map(x=>x["Species"].replace("Lys_pep_need_to_confirm_","")))).filter(x=>x!= "").filter(x=>x!=undefined)
                    let species_list = Array.from(new Set(data.map(x=>x["ID"].replace("Lys_pep_need_to_confirm_","")))).filter(x=>x!= "").filter(x=>x!=undefined)
                    let maxsymlen = calcuSymbollen(geneSymbol_list) + 5 //这个可以留
                    let maxproteinlen = calcuSymbollen(proteinid_list) + 5


                    // let specieslist = Array.from(new Set(data.map(x=>x["Species"].replace("Lys_pep_need_to_confirm_","")))).filter(x=>x!= "")
                    // let specieslist = Array.from(new Set(data.map(x=>x["ID"].replace("Lys_pep_need_to_confirm_","")))).filter(x=>x!= "")
                    
                    this.data.proteinid_list = proteinid_list
                    this.data.species_list = species_list
                    this.data.maxsymlen = maxsymlen
                    this.data.maxproteinlen = maxproteinlen
                    this.data.geneSymbol_list = geneSymbol_list
                    this.data.textSymbol_list = this.data.geneSymbol_list.map(x=>x.replace("gene_symbol:",""))
                    this.data.chosensymlist = this.data.geneSymbol_list.slice(0,5)
                    this.data.co2servermap = co2servermap

                    let colorMap = {}
                    species_list.forEach((ele,index) => {
                        const color = Oviz.color.Color.hsl((index%6)*60, 60+Math.floor((index/6))*10, 60+Math.floor((index/6))*10)
                        colorMap[ele] = color.string
                    });
                    this.data.colorMap = colorMap

                    this.data.showproteincolorMap = false
                    this.data.proteincolorMap = {}
                    this.data.showspeciescolorMap =false
                    this.data.speciescolorMap = {}
                    this.data.showspeciescolorMap =false
                    this.data.partscolorMap = {}
                    
                    this.data.showTitlelist = []
                    this.data.titlesep = getMaxlength(proteinid_list) + 10
                    this.data.geneMaxlen = getMaxlength(geneSymbol_list)


                    this.data.chosenGenesymbol = ""
                    this.data.selectgeneSymbol_list = {select:this.data.treefilslist.map(d=>d.replace("filter","").replace(".pep.aln.newick",""))} //这里和editor里改有关
                    
                    
                    this.data.geneBlockcolormap = {}
                    // debug 2023-03-31
                    let chromColorMap = {}
                    this.data.species_list.forEach((ele,index) => {
                        let chromosome = ele.split("-")[0]
                        if(chromColorMap[chromosome] == undefined){
                            const color = Oviz.color.Color.hsl((index%6)*60, 60+Math.floor((index/6))*10, 60+Math.floor((index/6))*10)
                            chromColorMap[chromosome] = color.string
                        }
                    })
                    //
                    this.data.geneSymbol_list.forEach((element,index) => {
                        // const color = Oviz.color.Color.hsl((index%6)*60, 60+Math.floor((index/6))*10, 60+Math.floor((index/6))*10)
                        // this.data.geneBlockcolormap[element] = color.string
                        // debug 2023-03-31
                        let chromosome = this.data.species_list[index].split("-")[0]
                        this.data.geneBlockcolormap[element] = chromColorMap[chromosome]
                        // debug 2023-03-31
                    });

                    // console.log("this.data.geneSymbol_list:",this.data.geneSymbol_list)
                    // console.log("this.data.geneBlockcolormap:",this.data.geneBlockcolormap)

                                                                       

        
                    // seq info data 和 tree 进行比较 找到交集 今儿确认叶子的数量
                    let infoData = {}
                    data.forEach(element => {
                        // infoData[element["Protein_ID"]] = element 
                        infoData[element["Best_Peptide"]] = element 
                    });
                    this.data.infoData = infoData

                    let childnodeslist = getList([this.data.genetree],this.data.infoData) //获得子节点的list
                    this.data.childnodeslist = childnodeslist

                    this.data.Leaveslist = []
                    getLeaveslist(this.data.genetree,this.data.Leaveslist) //用于计算扇形角度
                    //

                    // this.data.treeheight = 50 * 12 //这里需要改
                    this.data.treeheight = this.data.Leaveslist.length * 12

                    this.data.Leaveslist.length < 20? this.data.radicalR = this.data.Leaveslist.length * 20:this.data.radicalR = this.data.Leaveslist.length * 2
                    
                    this.data.treeColormap = {}
                    this.data.treeColormap["Potein Id Color"] = "#274851"
                    this.data.treeColormap["Node Color"] = "#809cb1"
                    this.data.treeColormap["Link Color"] = "black"

                }
            },
        },
        setup() {
            this.defineGradient("kg", "horizontal", ["#dbdbdb", "blue"]);
            Object.entries(this.data.geneBlockcolormap).forEach((item,index)=>{
                this.defineGradient(item[0], "horizontal", ["white", item[1]]);
            })   
            this.data.geneSymbol_list.forEach((element,index) => {
                const color = Oviz.color.Color.hsl((index%6)*60, 60+Math.floor((index/6))*10, 60+Math.floor((index/6))*10)
                this.defineGradient(element+"", "horizontal", ["white", color.string]);
                this.defineGradient(element+"Left", "horizontal", [color.string,"white"]);
            });
            this.size.width = this.data.treewidth + 1100
            this.size.height = this.data.treeheight + 200

            

            // registerEditorConfig(editorConfig(this), editorRef);
            registerEditorConfig(editorConfig(this), "getVue", "#task-output", editorRef);
            
        },
    });

    return visualizer;
}

export function registerAnalysisGeneTree() {
    register(MODULE_NAME, init);
}

register(MODULE_NAME, init);

export function getLeaveslist(node,arr){
    //For Oviz under the tree recursively find the list of leaf nodes
    node.children.forEach((item,index) => {
        Object.keys(item).indexOf("children") != -1? getLeaveslist(item,arr):arr.push(item.name)
    });
}

export function calcuSymbollen(symbollist){
    let lenlist = []
    symbollist.forEach(element => {
        lenlist.push(TextSize.measuredTextSize(element.replace("gene_symbol:",""),10).width)
    });
    return Math.max(...lenlist)
}

export function getColor(index){
    let colorlist = ["#031219","#045f72","#0a9496","#90d3c1","#ecd8a6","#ed9d00","#cc6601","#bb3f01","#ae1f11","#9f2027",
                    "#878576","#305659","#32202c","#9d382e","#a34f37","#a86344","#af9d85"]
    if(index>colorlist.length-1){
        index = index % colorlist.length
    }
    return colorlist[index]
}

export function convetTSV(oridata){
    let convertdata = {data:[],column:[]}
    let columname = oridata.split("\n")[0].split("\t")
    oridata.split("\n").forEach((element,index) => {
        let singlearray = {}
        element.split("\t").forEach((eachchr,chrindex) => {
            singlearray[columname[chrindex]] = eachchr
        });
        convertdata.data.push(singlearray)
    });
    convertdata.column = columname
    convertdata.data = convertdata.data.splice(1)
    return convertdata
}

 
let temp = [];
export function getList(tree,infoData) {
  for (let item in tree) {
    if (tree[item].children) {
      getList(tree[item].children,infoData)
    } else {
        // tree[item]["genesymbol"] = infoData[tree[item].name]["Gene_symbol"]
        // console.log(">>>>>>>DEBUG", infoData[tree[item].name])
        tree[item]["genesymbol"] = infoData[tree[item].name]["Gene"]
        temp.push(tree[item])
    }
  }
  return temp;
}
 
export const groupBy = (list, fn) => {
    const groups = {};
    list.forEach(function (o) {
        const group = JSON.stringify(fn(o));
        groups[group] = groups[group] || [];
        groups[group].push(o);
    });
    return groups;
}

export function getMaxlength(arr){
    let max = []
    arr.forEach(element => {
        max.push(TextSize.measuredTextSize(element,10).width)
    });
    return Math.max.apply(null,max)
}