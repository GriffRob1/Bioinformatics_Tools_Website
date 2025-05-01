import React, {useEffect} from 'react';
import {useState} from "react";
import {useNavigate} from "react-router-dom";

export default function ToolDescription({id, URL, imagePath, toolTitle, textDescription, isFavoritable=true}) {

    let favoritesList = localStorage.getItem('favorites')
    const [isFavorited, setIsFavorited] = useState((favoritesList.indexOf(id) >= 0));

    //adds favorited tools to localStorage and toggles isFavorited
    const onClickHeart = (e) => {
        e.stopPropagation();

        let newFavoriteArray = JSON.parse(localStorage.getItem('favorites'));
        if (isFavorited) {
            let index = newFavoriteArray.indexOf(id);
            newFavoriteArray.splice(index, 1);
        }
        else {
            newFavoriteArray.push(id)
        }
        localStorage.setItem('favorites', JSON.stringify(newFavoriteArray));

        setIsFavorited(!isFavorited);
    }

    const heartIconLink = isFavorited ? '/images/filled_in_heart.png' : '/images/empty_heart.png'

    const navigate = useNavigate();
    const navigateToToolPage = () => {
        navigate(URL);
    }

    return (
      <div className={'container tool-description'} onClick={navigateToToolPage}>
          {isFavoritable ? <img className={'heart-icon'} src={heartIconLink} alt={'heart'} onClick={onClickHeart}/> : ''}
          <img className={'tool-reference-image'} src={imagePath} alt={'description'}/>
          <h3 className={'tool-title'}>{toolTitle}</h3>
          <p className={'tool-text-description'}>{textDescription}</p>
      </div>
    );
}