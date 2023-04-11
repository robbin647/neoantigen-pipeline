import { CircosContentLine } from "crux/src/element/circos";
import { EditorDef } from "utils/editor";
// import { mkchildrentree } from "./data";
import { copyObject } from "utils/object";
import { getLeaveslist, getList } from "./index";

function run(v) {
    v.forceRedraw = true;
    v.run();
}

export const editorRef: any = {};

export const conf = {
    min: 0.00,
    treeDepth: 0,
    level: 7,
    nameOption: 2,
    isRadical: false,
    displayCircularTree: false,
    distinctNodeOnly: false,
    depthSelectOption: [
        { value: "0", text: "Root"},
        { value: "1", text: "Kingdom"},
        { value: "2", text: "Phylum"},
        { value: "3", text: "Gene Tree"},
        { value: "4", text: "Order"},
    ],
} as any;



export const generateTestConfig = (v): any => ({
    id: "plot_st",
    title: "General Setting",
    layout: "tabs",
    tabs:[
        {
            id:"g-common",
            name:"Size Setting",
            view:{
                type: "list",
                items:[
                    {
                        title: "Tree Width",
                        type: "input",
                        format: "int",
                        value: {
                            current: v.data.treewidth,
                            callback(x) {
                                v.data.treewidth = parseFloat(x);
                                v.size.width = v.data.treewidth + 1100
                                v.size.height = v.data.treeheight + 200
                                v.data._sizeUpdated = true;
                                run(v);
                            },
                        },
                    },
                    {
                        title: "Tree Height",
                        type: "input",
                        format: "int",
                        value: {
                            current: v.data.treeheight,
                            callback(x) {
                                v.data.treeheight = parseFloat(x);
                                v.size.width = v.data.treewidth + 1100
                                v.size.height = v.data.treeheight + 200
                                v.data._sizeUpdated = true;
                                run(v);
                            },
                        },
                    },
                ]

            }
        },
        
    ]
})

export const generateColorConfig = (v): any => ({
    id: "plot_st2",
    title: "Color Setting",
    layout: "tabs",
    tabs:[
        {
            id:"g-common",
            name:"Tree",
            view:{
                type: "list",
                items:[
                    {
                        title: "Apply ID Color Pattern",
                        type: "checkbox",
                        Option:v.data.specisColor,
                        value:{
                            current: v.data.specisColor,
                            callback(d){
                                v.data.specisColor = d
                                v.root._sizeUpdated = true;
                                v.forceRedraw = true;
                                run(v); 
                            }
                        }
                    },
                    {
                        type: "vue",
                        // title: "Method Color",
                        component: "color-picker",
                        data: {
                            title: "Node & Path Color",
                            scheme: copyObject(v.data.treeColormap),
                            id: "pwcolor",
                            callback(colors) {
                                v.data.treeColormap = colors;
                                v.data._sizeUpdated = true;
                                v.forceRedraw = true;
                                run(v);
                            },
                        },
                    },
                    {
                        type: "vue",
                        // title: "Method Color",
                        component: "color-picker",
                        data: {
                            title: "Leaves Color",
                            scheme: copyObject(v.data.colorMap),
                            id: "pwcolor",
                            callback(colors) {
                                v.data.colorMap = colors;
                                v.data._sizeUpdated = true;
                                v.forceRedraw = true;
                                run(v);
                            },
                        },
                    }
                ]

            }
        },
        {
            id:"g-commonklm",
            name:"Info Block",
            view:{
                type: "list",
                items:[
                    {
                        title: "Apply Color Pattern",
                        type: "checkbox",
                        Option:v.data.genecolorpattern,
                        value:{
                            current: v.data.genecolorpattern,
                            callback(d){
                                v.data.genecolorpattern = d
                                v.root._sizeUpdated = true;
                                v.forceRedraw = true;
                                run(v); 
                            }
                        }
                    },
                    // {
                    //     title: "Fill Fully",
                    //     type: "checkbox",
                    //     Option:v.data.geneBlockcolorfull,
                    //     value:{
                    //         current: v.data.geneBlockcolorfull,
                    //         callback(d){
                    //             v.data.geneBlockcolorfull = d
                    //             v.root._sizeUpdated = true;
                    //             v.forceRedraw = true;
                    //             run(v); 
                    //         }
                    //     }
                    // },
                    {
                        title: "Fill Opacity of Info Block",
                        type: "input",
                        value: {
                            current: v.data.fillOpacitynum,
                            callback(d) {
                                v.data.fillOpacitynum = parseFloat(d);
                                v.data._sizeUpdated = true;
                                v.forceRedraw = true;
                                run(v);
                            },
                        },
                    },
                    {
                        type: "vue",
                        component: "color-picker",
                        data: {
                            title: "Gene Block Color",
                            scheme: copyObject(v.data.geneBlockcolormap),
                            id: "pwcolor",
                            callback(colors) {
                                v.data.geneBlockcolormap = colors;
                                v.data._sizeUpdated = true;
                                v.forceRedraw = true;
                                run(v);
                            },
                        },
                    }
                ]

            }
        },
        
    ]
})

export const generateCusConfig = (v):any => ({
    id: "plot-stb",
    title: "Customized Setting",
    layout: "tabs",
    tabs: [
        // [DEBUG 2023-03-31] The function that this part is to do is meaningless now, so I comment it out.
        // {
        //     id:"sng"+"csty-common",
        //     name:"Setting",
        //     view:{
        //         type:"list",
        //         items:[
        //             // {
        //             //     title: "Protein Id",
        //             //     type: "checkbox",
        //             //     Option:v.data.showprotein,
        //             //     value:{
        //             //         callback(d){
        //             //             v.data.showprotein = d
        //             //             v.data.showprotein? v.data.showTitlelist.push("Protein id"):v.data.showTitlelist.splice(v.data.showTitlelist.indexOf("Protein id"),1)
        //             //             v.root._sizeUpdated = true;
        //             //             v.forceRedraw = true;
        //             //             run(v); 
        //             //         }
        //             //     }
        //             // },
        //             {
        //                 title: "Leaf Node Separator",
        //                 type: "input",
        //                 value: {
        //                     current: v.data.leafSeparator,
        //                     callback(d) {
        //                         v.data.leafSeparator = d+"";
        //                         v.data._sizeUpdated = true;
        //                         v.forceRedraw = true;
        //                         run(v);
        //                     },
        //                 },
        //             },
        //         ]
        //     }
        // },
        // [DEBUG 2023-03-31] This part is not working. Omit it for now.
        // {
        //     id:"sngnn"+"csty-common",
        //     name:"Gene setting",
        //     view:{
        //         type:"list",
        //         items:[
        //             {
        //                 title: "Choose Genesymbol",
        //                 type: "select",
        //                 ref: "depthSelect",
        //                 options: v.data.selectgeneSymbol_list.select,
        //                 bind: {
        //                     object: v.data.selectgeneSymbol_list,
        //                     path: "select",
        //                     callback(d) {
        //                         // editorRef.oscolor.config.data.scheme = v.data.OSPlotdata.colormap
        //                         v.data.chosenGenesymbol = d
        //                         console.log("v.data.chosenGenesymbol:",v.data.chosenGenesymbol)
        //                         console.log("v.data.chosenGenesymbol:",v.data.chosenGenesymbol)
        //                         mkchildrentree(v)
        //                         v.run();
        //                     },
        //                 },
        //             },
        //         ]
        //     }
        // },
        {
            id:"sng23"+"csty-common",
            name:"Formal Setting",
            view:{
                type:"list",
                items:[
                    {
                        title: "Choose Horizontal or Radical",
                        type: "select",
                        options: ["Horizontal","Radical"],
                        value:{
                            callback(d) {
                                v.data.horizonalPart = d == "Horizontal"? true:false
                                v.data.radicalPart = d == "Radical"? true:false
                                v.root._sizeUpdated = true;
                                v.forceRedraw = true;
                                run(v);
                            }
                        }
                    },
                    // {
                    //     title: "Horizontal",
                    //     type: "checkbox",
                    //     Option:v.data.horizonalPart,
                    //     value:{
                    //         callback(d){
                    //             v.data.horizonalPart = d
                    //             v.data.radicalPart = !d
                    //             v.root._sizeUpdated = true;
                    //             v.forceRedraw = true;
                    //             run(v); 
                    //         }
                    //     }
                    // },
                    // {
                    //     title: "Radical",
                    //     type: "checkbox",
                    //     Option:v.data.radicalPart,
                    //     value:{
                    //         callback(d){
                    //             v.data.radicalPart = d
                    //             v.data.horizonalPart = !d
                    //             console.log("radicalR:",v.data.radicalR)
                    //             v.size.width = 2*(v.data.radicalR*1.5 + 230 + v.data.geneMaxlen + v.data.maxproteinlen + 10)
                    //             v.size.height = 2*(v.data.radicalR + 230 + v.data.geneMaxlen + v.data.maxproteinlen + 10)
                    //             v.root._sizeUpdated = true;
                    //             v.forceRedraw = true;
                    //             run(v); 
                    //         }
                    //     }
                    // },
                ]
            }
        },
    ]
})


export function mkchildrentree(v){
    // v.data.genetree = v.data.allDatas[v.data.oritreefilslist.indexOf(v.data.chosenGenesymbol)]
    console.log("before_v.data.genetree:",v.data.genetree)
    console.log("v.data.allDatas:",v.data.allDatas)
    console.log("oritreefilslist:",v.data.oritreefilslist)
    // console.log("index:",v.data.oritreefilslist.indexOf(v.data.chosenGenesymbol))
    v.data.oritreefilslist.forEach((item,index) => {
        item.indexOf(v.data.chosenGenesymbol)!= -1? v.data.genetree = v.data.allDatas[index]:null
    });
    console.log("v.data.genetree:",v.data.genetree)

    let testnodesList = []
    childNodesDeep(v.data.genetree.children,testnodesList)
    v.data.treeheight = testnodesList.length * 12  
    v.size.width = v.data.treewidth + 1100
    v.size.height = v.data.treeheight + 200
    v.data.genecolorpattern = false
    v.data.Leaveslist = []
    getLeaveslist(v.data.genetree,v.data.Leaveslist) //用于计算扇形角度
    v.data.Leaveslist.length < 20? v.data.radicalR = v.data.Leaveslist.length * 20:v.data.radicalR = v.data.Leaveslist.length * 2
                    
}

//获得最后一层的节点
export function childNodesDeep(nodes,arr){
    nodes.forEach((ele,index) => {
        arr.push(ele.name);
        if(ele.children){
            childNodesDeep(ele.children,arr);
        }
    })
}


export function lookForAllId(infoData,data = [], arr = [],ob = {}) {
    for (let item of data) {
        if(infoData[item.name]) arr.push(infoData[item.name].gene_symbol)
        if (item.children && item.children.length) lookForAllId(infoData,item.children, arr,ob)
    }
    return arr
}





export function editorConfig(v): EditorDef {
    return {
        sections: [
            generateTestConfig(v),
            generateColorConfig(v),
            generateCusConfig(v),
        ]
    }
}