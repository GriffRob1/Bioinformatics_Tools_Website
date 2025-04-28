import React from 'react';

export default function ToolsSearchBar() {
    return (
        <form className={'container tools-search-bar'} action={'/tools'} method={'GET'}>
            <fieldset className={'container filter-section'}>
                <legend>Filter Tools:</legend>

                <label className={'tool-name-input'}>
                    Name of Tool:
                    <input type={'text'} name={'name'} id={'name'} />
                </label>
                <label className={'category-input'}>
                    Tool Category:
                    <select name={'category'} id={'category'}>
                        <option disabled selected>Enter Category</option>
                        <option>Motif Finding</option>
                        <option>Genome Sequencing</option>
                        <option>Pairwise Alignment</option>
                    </select>
                </label>

                <label className={'show-favorites-input'}>
                    Show only favorites:
                    <input type={'checkbox'} name={'show-favorites'} id={'show-favorites'}/>
                </label>
            </fieldset>
            <fieldset className={'container sort-section'}>
                <legend>Sort By:</legend>

                <div className={'container sorting-radio-buttons'}>
                    <label className={'sort-by-most-popular-input'}>
                        Most Popular
                        <input type={'radio'} name={'sort-by'} value={'sort-by-most-popular'}/>
                    </label>

                    <label className={'sort-by-alphabetical-input'}>
                        Alphabetical
                        <input type={'radio'} name={'sort-by'} value={'sort-by-alphabetical'}/>
                    </label>

                    <label className={'sort-by-newest-input'}>
                        Newest
                        <input type={'radio'} name={'sort-by'} value={'sort-by-newest'}/>
                    </label>
                </div>
            </fieldset>
            <button type={'submit'} className={'blue-button search-button'}>Search Tools</button>
        </form>
    )
}