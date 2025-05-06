import React from 'react';
import {useState} from "react";
import ToolsSearchBar from "../components/ToolsSearchBar";
import ToolsSearchContainer from "../components/ToolsSearchContainer";

export default function Tools({toolsList}) {
    const [inputName, setInputName] = useState('');
    const [inputCategory, setInputCategory] = useState('');
    const [inputShowOnlyFavorites, setInputShowOnlyFavorites] = useState(false);
    const [inputSortBy, setInputSortBy] = useState('sort-by-most-popular');

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
            <ToolsSearchContainer toolsList={toolsList}
                                  inputName={inputName}
                                  inputCategory={inputCategory}
                                  inputSortBy={inputSortBy}
                                  inputShowOnlyFavorites={inputShowOnlyFavorites}/>
        </section>
    )
}