import React from 'react';
import {useState} from "react";
import ToolsSearchBar from "../components/ToolsSearchBar";
import ToolsSearchContainer from "../components/ToolsSearchContainer";

export default function Tools() {

    let toolsSampleData = `[
        {
            "id": 0,
            "popularity": 9,
            "URL": "/genome-assembler",
            "imagePath": "/images/favicon.ico",
            "toolTitle": "Genome Assembler",
            "textDescription": "Assemble a genome using short nucleotide reads",
            "dateAdded": "2025-01-17T03:24:00",
            "category": "Genome Sequencing Tools"
        },

        {
            "id": 1,
            "popularity": 5,
            "URL": "/motif-finder",
            "imagePath": "/images/logo192.png",
            "toolTitle": "Motif Finder",
            "textDescription": "Find a common pattern in a list of reads",
            "dateAdded": "2025-01-18T03:24:00",
            "category": "Motif Finding Tools"
        },

        {
            "id": 2,
            "popularity": 3,
            "URL": "/sequence-alignment",
            "imagePath": "/images/logo512.png",
            "toolTitle": "Sequence Alignment",
            "textDescription": "Align 2 or more reads using an alignment algorithm",
            "dateAdded": "2025-01-19T03:24:00",
            "category": "Sequence Alignment Tools"
        },
        
        {
            "id": 3,
            "popularity": 9,
            "URL": "/genome-assembler",
            "imagePath": "/images/favicon.ico",
            "toolTitle": "Genome Assembler",
            "textDescription": "Assemble a genome using short nucleotide reads",
            "dateAdded": "2025-01-17T03:24:00",
            "category": "Genome Sequencing Tools"
        },

        {
            "id": 4,
            "popularity": 5,
            "URL": "/motif-finder",
            "imagePath": "/images/logo192.png",
            "toolTitle": "Motif Finder",
            "textDescription": "Find a common pattern in a list of reads",
            "dateAdded": "2025-01-18T03:24:00",
            "category": "Motif Finding Tools"
        },

        {
            "id": 5,
            "popularity": 3,
            "URL": "/sequence-alignment",
            "imagePath": "/images/logo512.png",
            "toolTitle": "Sequence Alignment",
            "textDescription": "Align 2 or more reads using an alignment algorithm",
            "dateAdded": "2025-01-19T03:24:00",
            "category": "Sequence Alignment Tools"
        },
        
        {
            "id": 6,
            "popularity": 9,
            "URL": "/genome-assembler",
            "imagePath": "/images/favicon.ico",
            "toolTitle": "Genome Assembler",
            "textDescription": "Assemble a genome using short nucleotide reads",
            "dateAdded": "2025-01-17T03:24:00",
            "category": "Genome Sequencing Tools"
        },

        {
            "id": 7,
            "popularity": 5,
            "URL": "/motif-finder",
            "imagePath": "/images/logo192.png",
            "toolTitle": "Motif Finder",
            "textDescription": "Find a common pattern in a list of reads",
            "dateAdded": "2025-01-18T03:24:00",
            "category": "Motif Finding Tools"
        },

        {
            "id": 8,
            "popularity": 3,
            "URL": "/sequence-alignment",
            "imagePath": "/images/logo512.png",
            "toolTitle": "Sequence Alignment",
            "textDescription": "Align 2 or more reads using an alignment algorithm",
            "dateAdded": "2025-01-19T03:24:00",
            "category": "Sequence Alignment Tools"
        }
    ]`

    const [inputName, setInputName] = useState('');
    const [inputCategory, setInputCategory] = useState('');
    const [inputSortBy, setInputSortBy] = useState('sort-by-most-popular');

    return (
        <section className={'container tools-page'}>
            <h1 className={'search-tools-title'}>Search Our Tools</h1>
            <ToolsSearchBar toolsSampleData={toolsSampleData}
                            inputName={inputName}
                            inputCategory={inputCategory}
                            inputSortBy={inputSortBy}
                            onInputNameChange={setInputName}
                            onInputCategoryChange={setInputCategory}
                            onInputSortByChange={setInputSortBy}/>
            <ToolsSearchContainer toolsSampleData={toolsSampleData}
                                  inputName={inputName}
                                  inputCategory={inputCategory}
                                  inputSortBy={inputSortBy}/>
        </section>
    )
}