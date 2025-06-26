import '../styles/Tools.css';
import React from 'react';
import {useState} from "react";
import ToolsSearchBar from "../components/ToolsSearchBar";
import ToolsSearchContainer from "../components/ToolsSearchContainer";
import BlueButton from "../components/BlueButton";

export default function Tools({toolsList, setToolsList}) {
    const [inputName, setInputName] = useState('');
    const [inputCategory, setInputCategory] = useState('');
    const [inputShowOnlyFavorites, setInputShowOnlyFavorites] = useState(false);
    const [inputSortBy, setInputSortBy] = useState('sort-by-newest');

    return (
        <section className={'container tools-page'}>
            <h1 className={'search-tools-title'}>Search Our Tools</h1>
            <ToolsSearchBar inputName={inputName}
                            inputCategory={inputCategory}
                            inputShowOnlyFavorites={inputShowOnlyFavorites}
                            inputSortBy={inputSortBy}
                            onInputNameChange={setInputName}
                            onInputCategoryChange={setInputCategory}
                            onInputShowOnlyFavoritesChange={setInputShowOnlyFavorites}
                            onInputSortByChange={setInputSortBy}/>
            <ToolsSearchContainer inputName={inputName}
                                  inputCategory={inputCategory}
                                  inputSortBy={inputSortBy}
                                  inputShowOnlyFavorites={inputShowOnlyFavorites}
                                  toolsList={toolsList}
                                  setToolsList={setToolsList}/>
            <BlueButton onClick={() => window.scrollTo(0,0)} buttonClass={'scroll-to-top-button'}>Scroll to top</BlueButton>
        </section>
    )
}