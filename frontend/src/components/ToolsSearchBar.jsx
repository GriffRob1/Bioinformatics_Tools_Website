import '../styles/ToolsSearchBar.css';
import React from 'react';

export default function ToolsSearchBar({inputName,
                                        inputCategory,
                                        inputShowOnlyFavorites,
                                        inputSortBy,
                                        onInputNameChange,
                                        onInputCategoryChange,
                                        onInputShowOnlyFavoritesChange,
                                        onInputSortByChange}) {
    return (
        <form className={'container tools-search-bar'} onSubmit={(e) => e.preventDefault()}>
            <fieldset className={'container filter-section'}>
                <legend>Filter Tools:</legend>

                <label className={'tool-name-input'}>
                    Name of Tool:
                    <input type={'text'}
                           name={'name'}
                           id={'name'}
                           value={inputName}
                           onChange={(e) => onInputNameChange(e.target.value)}/>
                </label>
                <label className={'category-input'}>
                    Tool Category:
                    <select
                        name={'category'}
                        id={'category'}
                        value={inputCategory}
                        onChange={(e) => onInputCategoryChange(e.target.value)}
                    >
                        <option value={''}>No Category Selected</option>
                        <option value={'Frequent Pattern Finding'}>Frequent Pattern Finding</option>
                        <option value={'Motif Finding'}>Motif Finding</option>
                        <option value={'Genome Assembly'}>Genome Assembly</option>
                        <option value={'Cyclopeptide Sequencing'}>Cyclopeptide Sequencing</option>
                    </select>
                </label>

                <label className={'show-favorites-input'}>
                    Show only favorites:
                    <input type={'checkbox'}
                           name={'show-favorites'}
                           id={'show-favorites'}
                           onChange={(e) => onInputShowOnlyFavoritesChange(e.target.checked)}/>
                </label>
            </fieldset>
            <fieldset className={'container sort-section'}>
                <legend>Sort By:</legend>

                <div className={'container sorting-radio-buttons'}>

                    <label className={'sort-by-most-popular-input'}>
                        Most Popular
                        <input type={'radio'}
                               name={'sort-by'}
                               value={'sort-by-most-popular'}
                               onChange={(e) => onInputSortByChange(e.target.value)}/>
                    </label>

                    <label className={'sort-by-alphabetical-input'}>
                        Alphabetical
                        <input type={'radio'}
                               name={'sort-by'}
                               value={'sort-by-alphabetical'}
                               onChange={(e) => onInputSortByChange(e.target.value)}/>
                    </label>

                    <label className={'sort-by-newest-input'}>
                        Newest
                        <input type={'radio'}
                               name={'sort-by'}
                               value={'sort-by-newest'}
                               checked={inputSortBy === 'sort-by-newest'}
                               onChange={(e) => onInputSortByChange(e.target.value)}/>
                    </label>
                </div>
            </fieldset>
        </form>
    )
}